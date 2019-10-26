const [N, S, E, W, NE, NW, SE, SW] = [1, 2, 3, 4, 5, 6, 7, 8];
const opp = [0, S, N, W, E, SW, SE, NW, NE];
const weights = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
const cxs = [0, 0, 0, 1, -1, 1, -1, 1, -1];
const cys = [0, 1, -1, 0, 0, 1, 1, -1, -1];

function createVectors(size: number) {
  const buffer = new ArrayBuffer(4*9*size);
  return Array.from({ length: size }, (_, i) => {
    const v = new Float32Array(buffer, 4*9*i, 9);
    v.set(weights);
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
    const size = xdim * ydim;

    this.vectors = createVectors(size);
    this.swap = createVectors(size);
  
    this.rho = new Float32Array(size);
    this.rho.fill(1);
  }

  // Collide particles within each cell (here's the physics!):
  public collide(viscosity: number) {
    const { xdim, ydim, rho, vectors } = this;
    const omega = 1 / (3*viscosity + 0.5); // reciprocal of relaxation time
    const invomega = 1 - omega;
    const max = xdim * ydim;
    for (let i=0; i<max; i++) {
      // vectors for this lattice location
      const vecs = vectors[i];

      // macroscopic density
      const newrho = vecs.reduce((a, v) => a + v);
      rho[i] = newrho;
      
      // macroscopic velocity components
      let ux = 0, uy = 0;
      for (let j = 1; j < 9; j++) {
        const v = vecs[j];
        ux += cxs[j]*v
        uy += cys[j]*v;
      }
      ux/=newrho;
      uy/=newrho;

      /*
      f_j^eq = w_j rho (1 + (e_j|u) / c_s^2 + (e_j|u)^2/(2 c_s^4) - u^2/(2 c_s^2))
      f`_j = omega f_j^eq + (1 - omega) f_j 

      c_s, the lattice speed of sound, is 1/sqrt(3)
      Thus, 1/c_s^2 = 3
            1/(2 c_s^4) = 4.5
            1/(2 c_s^2) = 1.5
      */

      const u2 =  1 - 1.5 * (ux * ux + uy * uy);
      for (let j = 0; j < 9; j++) {
        const dir = cxs[j]*ux + cys[j]*uy;
        const eq = weights[j] * newrho * (u2 + 3 * dir + 4.5 * dir * dir);
        vecs[j] = omega * eq + invomega*vecs[j];
      }
    }

    //this.swap = vectors;
    //this.vectors = swap;
  }

  // Move particles along their directions of motion:
  public stream(barriers: boolean[]) {
    const { xdim, ydim, vectors, swap } = this;

    const size = xdim*ydim;
    for (let i = 1; i < size; i++) {
      swap[i].set(vectors[i]);
    }

    const index = (x: number, y: number) => (x%xdim)+(y%ydim)*xdim;

    for (let y=1; y<ydim-1; y++) {
      for (let x=1; x<xdim-1; x++) {
        const i = index(x, y);
        const v = vectors[i];
        for (let j=1;j<9;j++) {
          v[j] = swap[index(x-cxs[j], y-cys[j])][j]; // move particles
        }
      }
    }
    for (let y=0; y<ydim; y++) { // Now handle bounce-back from barriers
      for (let x=0; x<xdim; x++) {
        const i = index(x, y);
        if (barriers[i]) {
          const s = swap[i];
          for (let j=1;j<9;j++) {
            vectors[index(x+cxs[j], y+cys[j])][j] = s[opp[j]]; // move particles
          }
          // Force on the barrier:
          // barrierFx += v[E][i] + v[NE][i] + v[SE][i] - v[W][i] - v[NW][i] - v[SW][i];
          // barrierFy += v[N][i] + v[NE][i] + v[NW][i] - v[S][i] - v[SE][i] - v[SW][i];
        }
      }
    }
  }

  // Set all densities in a cell to their equilibrium values for a given velocity and density:
  public setEquilibrium(x: number, y: number, ux: number, uy: number, rho: number) {
    const { xdim, vectors } = this;
    const i = x + y*xdim;
    this.rho[i] = rho;
    const v = vectors[i];

    const u2 =  1.5 * (ux * ux + uy * uy);
    for (let j = 0; j < 9; j++) {
      const dir = cxs[j]*ux + cys[j]*uy;
      v[j] =  weights[j] * rho * (1 + 3 * dir + 4.5 * dir * dir - u2);
    }
  }
}