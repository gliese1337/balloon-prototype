const [N, S, E, W, NE, NW, SE, SW] = [1, 2, 3, 4, 5, 6, 7, 8];
const opp = [0, S, N, W, E, SW, SE, NW, NE];
const weights = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
const cxs = [0, 0, 0, 1, -1, 1, -1, 1, -1];
const cys = [0, 1, -1, 0, 0, 1, 1, -1, -1];

function createVectors(size: number) {
  const wordsize = 9*size;
  const buffer = new Float32Array(wordsize);
  for (let i = 0; i < wordsize; i+=9) {
    buffer.set(weights, i);
  }
  return buffer;
}

export default class LatticeBoltzmann {        
  private streamed: Float32Array; // microscopic densities along each lattice direction
  private collided: Float32Array;
  public rho: Float32Array;    // macroscopic density; cached for rendering

  constructor(public xdim: number, public ydim: number) {
    const size = xdim * ydim;
    this.streamed = createVectors(size);
    this.collided = createVectors(size);
    this.rho = new Float32Array(size);
  }

  public collide(viscosity: number) {
    const { xdim, ydim, rho, collided, streamed } = this;
    const tau = 3*viscosity + 0.5; // relation timescale
    const omega = 1 / tau;
    const invomega = 1 - omega;
    const max = xdim * ydim;
    for (let i=0,i9=0; i<max; i++,i9+=9) {
      // macroscopic density
      let newrho = 0;
      for (let j=0;j<9;j++){
        newrho += streamed[i9+j];
      }
      rho[i] = newrho;
      
      // macroscopic velocity components
      let ux = 0, uy = 0;
      for (let j = 1; j < 9; j++) {
        const v = streamed[i9+j];
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
        collided[i9+j] = omega * eq + invomega * streamed[i9+j];
      }
    }
  }

  public stream(barriers: boolean[]) {
    const { xdim, ydim, collided, streamed } = this;
    const index = (x: number, y: number) => (x%xdim)+(y%ydim)*xdim;
    const cIndex = (x: number, y: number, s: -1|1, j: number) =>
                      9*(((x+s*cxs[j])%xdim)+((y+s*cys[j])%ydim)*xdim)+j;

    // Move particles along their directions of motion:
    for (let y=1; y<ydim-1; y++) {
      for (let x=1; x<xdim-1; x++) {
        const i = index(x, y);
        const i9 = i*9;
        for (let j=0;j<9;j++) {
          streamed[i9 + j] = collided[cIndex(x, y, -1, j)];
        }
      }
    }

    // Handle bounce-back from barriers
    for (let y=0; y<ydim; y++) {
      for (let x=0; x<xdim; x++) {
        const i = index(x, y);
        const i9 = i*9;
        if (barriers[i]) {
          for (let j=1;j<9;j++) {
            streamed[cIndex(x, y, 1, j)] = collided[i9 + opp[j]];
          }
        }
      }
    }
  }

  // Set all densities in a cell to their equilibrium values for a given velocity and density:
  public setEquilibrium(x: number, y: number, ux: number, uy: number, rho: number) {
    const { xdim, streamed } = this;
    const i = x + y*xdim;
    this.rho[i] = rho;

    const i9 = i*9;
    const u2 =  1 - 1.5 * (ux * ux + uy * uy);
    for (let j = 0; j < 9; j++) {
      const dir = cxs[j]*ux + cys[j]*uy;
      streamed[i9+j] =  weights[j] * rho * (u2 + 3 * dir + 4.5 * dir * dir);
    }
  }
}