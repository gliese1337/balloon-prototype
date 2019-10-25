// [Z, N, S, E, W, NE, NW, SE, SW]
const weights = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
const cxs = [0, 0, 0, 1, -1, 1, -1, 1, -1];
const cys = [0, 1, -1, 0, 0, 1, 1, -1, -1];

export default class LatticeBoltzmann {
  private vectors: Float32Array[];            // microscopic densities along each lattice direction
  public rho: Float32Array;           // macroscopic density

  constructor(public xdim: number, public ydim: number) {
    // Create the arrays of fluid particle densities, etc. (using 1D arrays for speed):
    // To index into these arrays, use x + y*xdim, traversing rows first and then columns.
    const size = 4*xdim * ydim;

    const buffer = new ArrayBuffer(40*size);
    const vectors = Array.from({ length: 9 }, (_, i) => {
      const v = new Float32Array(buffer, 4*i*size, size);
      v.fill(weights[i]);
      return v;
    });

    this.vectors = vectors;
  
    this.rho = new Float32Array(buffer, 36*size, size);
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
  }

  // Move particles along their directions of motion:
  public stream(barriers: boolean[]) {
    const { xdim, ydim, vectors } = this;
    const [, N, S, E, W, NE, NW, SE, SW] = vectors;

    const max = xdim * ydim;
    for (let y=xdim-1; y<max; y+=xdim) {
      // at right end, copy left-flowing densities from next column to the left
      W[y] = W[y-1];
      NW[y] = NW[y-1];
      SW[y] = SW[y-1];
    }

    for (let y=(ydim-2)*xdim; y>0; y-=xdim) {
      const xlimit = xdim + y - 1;
      for (let x=1+y; x<xlimit; x++) {            // first start in NW corner...
        N[x] = N[x-xdim];            // move the north-moving particles
        NW[x] = NW[x+1-xdim];        // and the northwest-moving particles
      }
      for (let x=xdim-2+y; x>y; x--) {            // now start in NE corner...
        E[x] = E[x-1];            // move the east-moving particles
        NE[x] = NE[x-1-xdim];        // and the northeast-moving particles
      }
    }
    const ylimit = (ydim-1)*xdim;
    for (let y=xdim; y<ylimit; y+=xdim) {
      for (let x=xdim-2+y; x>y; x--) {             // now start in SE corner...
        S[x] = S[x+xdim];            // move the south-moving particles
        SE[x] = SE[x-1+xdim];        // and the southeast-moving particles
      }
      const xlimit = xdim + y - 1;
      for (let x=1+y; x<xlimit; x++) { // now start in the SW corner...
        W[x] = W[x+1];            // move the west-moving particles
        SW[x] = SW[x+1+xdim];        // and the southwest-moving particles
      }
    }
    for (let y=xdim; y<ylimit; y+=xdim) {                // Now handle bounce-back from barriers
      for (let x=1+y; x<y+xdim-1; x++) {
        if (barriers[x]) {
            E[x+1] = W[x];
            W[x-1] = E[x];
            N[x+xdim] = S[x];
            S[x-xdim] = N[x];
            NE[x+1+xdim] = SW[x];
            NW[x-1+xdim] = SE[x];
            SE[x+1-xdim] = NW[x];
            SW[x-1-xdim] = NE[x];
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