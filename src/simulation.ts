// [Z, N, S, E, W, NE, NW, SE, SW]
const weights = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
const cxs = [0, 0, 0, 1, -1, 1, -1, 1, -1];
const cys = [0, 1, -1, 0, 0, 1, 1, -1, -1];

function createVectors(size: number) {
  const buffer = new ArrayBuffer(36*size);
  return Array.from({ length: 9 }, (_, i) => {
    const v = new Float32Array(buffer, 4*i*size, size);
    v.fill(weights[i]);
    return v;
  });
}

export default class LatticeBoltzmann {
  private vectors: Float32Array[];            // microscopic densities along each lattice direction
  private swap: Float32Array[];
  public rho: Float32Array;           // macroscopic density

  constructor(public xdim: number, public ydim: number) {
    // Create the arrays of fluid particle densities, etc. (using 1D arrays for speed):
    // To index into these arrays, use x + y*xdim, traversing rows first and then columns.
    const size = 4*xdim * ydim;

    this.vectors = createVectors(size);
    this.swap = createVectors(size);
  
    this.rho = new Float32Array(size);
    this.rho.fill(1);
  }

  // Collide particles within each cell (here's the physics!):
  public collide(viscosity: number) {
    const { xdim, ydim, rho, vectors } = this;
    const [Z, N, S, E, W, NE, NW, SE, SW] = vectors;
    const omega = 1 / (3*viscosity + 0.5); // reciprocal of relaxation time
    const invomega = 1 - omega;
    const max = xdim * ydim;
    for (let i=0; i<max; i++) {
      const newrho = Z[i]+N[i]+S[i]+E[i]+W[i]+NW[i]+NE[i]+SW[i]+SE[i];
      rho[i] = newrho;

      /*
      f_j^eq = w_j rho (1 + (e_j|u) / c_s^2 + (e_j|u)^2/(2 c_s^4) - u^2/(2 c_s^2))
      f`_j = omega f_j^eq + (1 - omega) f_j 

      c_s, the lattice speed of sound, is 1/sqrt(3)
      Thus, 1/c_s^2 = 3
            1/(2 c_s^4) = 4.5
            1/(2 c_s^2) = 1.5
      */

      // macroscopic velocity components
      const ux = (E[i]+NE[i]+SE[i]-W[i]-NW[i]-SW[i])/newrho;
      const uy = (N[i]+NE[i]+NW[i]-S[i]-SE[i]-SW[i])/newrho;
      const u2 =  1.5 * (ux * ux + uy * uy);
      for (let j = 0; j < 9; j++) {
        const v = vectors[j];
        const dir = cxs[j]*ux + cys[j]*uy;
        const eq = weights[j] * newrho * (1 + 3 * dir + 4.5 * dir * dir - u2);
        v[i] = omega * eq + invomega*v[i];
      }
    }

    //this.swap = vectors;
    //this.vectors = swap;
  }

  // Move particles along their directions of motion:
  public stream(barriers: boolean[]) {
    const { xdim, ydim, vectors, swap } = this;
    const [, N, S, E, W, NE, NW, SE, SW] = vectors;
    const [, sN, sS, sE, sW, sNE, sNW, sSE, sSW] = swap;

    sN.set(N); sS.set(S); sE.set(E); sW.set(W);
    sNE.set(NE); sNW.set(NW); sSE.set(SE); sSW.set(SW);

    const index = (x: number, y: number) => (x%xdim)+(y%ydim)*xdim;

    for (let y=1; y<ydim-1; y++) {
      for (let x=1; x<xdim-1; x++) {
        const i = index(x, y);
        N[i] = sN[index(x, y-1)];     // move the north-moving particles
        NW[i] = sNW[index(x+1, y-1)]; // move the northwest-moving particles
        E[i] = sE[index(x-1, y)];     // move the east-moving particles
        NE[i] = sNE[index(x-1, y-1)]; // move the northeast-moving particles
        S[i] = sS[index(x, y+1)];     // move the south-moving particles
        SE[i] = sSE[index(x-1, y+1)]; // move the southeast-moving particles
        W[i] = sW[index(x+1, y)];     // move the west-moving particles
        SW[i] = sSW[index(x+1, y+1)]; // move the southwest-moving particles
      }
    }
    for (let y=0; y<ydim; y++) { // Now handle bounce-back from barriers
      for (let x=0; x<xdim; x++) {
        const i = index(x, y);
        if (barriers[i]) {
            E[index(x+1, y)] = W[i];
            W[index(x-1, y)] = E[i];
            N[index(x, y+1)] = S[i];
            S[index(x, y-1)] = N[i];
            NE[index(x+1, y+1)] = SW[i];
            NW[index(x-1, y+1)] = SE[i];
            SE[index(x+1, y-1)] = NW[i];
            SW[index(x-1, y-1)] = NE[i];
            // Force on the barrier:
            // barrierFx += nE[index] + nNE[index] + nSE[index] - nW[index] - nNW[index] - nSW[index];
            // barrierFy += nN[index] + nNE[index] + nNW[index] - nS[index] - nSE[index] - nSW[index];
        }
      }
    }
  }

  // Set all densities in a cell to their equilibrium values for a given velocity and density:
  public setEquilibrium(x: number, y: number, ux: number, uy: number, rho: number) {
    const { xdim, vectors } = this;
    const i = x + y*xdim;
    this.rho[i] = rho;

    const u2 =  1.5 * (ux * ux + uy * uy);
    for (let j = 0; j < 9; j++) {
      const dir = cxs[j]*ux + cys[j]*uy;
      vectors[j][i] =  weights[j] * rho * (1 + 3 * dir + 4.5 * dir * dir - u2);
    }
  }
}