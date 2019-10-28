const weights = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
const cxs =     [  0,   0,   0,   1,  -1,    1,   -1,    1,   -1];
const cys =     [  0,   1,  -1,   0,   0,    1,   -1,   -1,    1];
const opp =     [  0,   2,   1,   4,   3,    6,    5,    8,    7];

const q = 9;

function createVectors(size: number) {
  const wordsize = q*size;
  const buffer = new Float32Array(wordsize);
  for (let i = 0; i < wordsize; i+=q) {
    buffer.set(weights, i);
  }
  return buffer;
}

export default class LatticeBoltzmann {        
  private streamed: Float32Array; // microscopic densities along each lattice direction
  private collided: Float32Array;
  private _rho: Float32Array;    // macroscopic density; cached for rendering

  constructor(public readonly xdim: number, public readonly ydim: number) {
    const size = xdim * ydim;
    this.streamed = createVectors(size);
    this.collided = createVectors(size);
    this._rho = new Float32Array(size);
  }

  public rho(x: number, y: number) {
    return this._rho[(x%this.xdim)+(y%this.ydim)*this.xdim];
  }

  public index(x: number, y: number) {
    return (x%this.xdim)+(y%this.ydim)*this.xdim;
  }

  public collide(viscosity: number) {
    const { xdim, ydim, _rho, collided, streamed } = this;
    const tau = 3*viscosity + 0.5; // relation timescale
    const omega = 1 / tau;
    const invomega = 1 - omega;
    const max = xdim * ydim;
    for (let i=0,iq=0; i<max; i++,iq+=q) {
      // macroscopic density
      let newrho = 0;
      for (let j=0;j<q;j++){
        newrho += streamed[iq+j];
      }
      _rho[i] = newrho;
      
      // macroscopic velocity components
      let ux = 0, uy = 0;
      for (let j = 1; j < q; j++) {
        const v = streamed[iq+j];
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
      for (let j = 0; j < q; j++) {
        const dir = cxs[j]*ux + cys[j]*uy;
        const eq = weights[j] * newrho * (u2 + 3 * dir + 4.5 * dir * dir);
        collided[iq+j] = omega * eq + invomega * streamed[iq+j];
      }
    }
  }

  public stream(barriers: boolean[]) {
    const { xdim, ydim, collided, streamed } = this;
    const index = (x: number, y: number) => (x%xdim)+(y%ydim)*xdim;
    const cIndex = (x: number, y: number, s: -1|1, j: number) =>
                      q*(((x+s*cxs[j])%xdim)+((y+s*cys[j])%ydim)*xdim)+j;

    // Move particles along their directions of motion:
    for (let y=1; y<ydim-1; y++) {
      for (let x=1; x<xdim-1; x++) {
        const i = index(x, y);
        const iq = i*q;
        for (let j=0;j<q;j++) {
          streamed[iq + j] = collided[cIndex(x, y, -1, j)];
        }
      }
    }

    // Handle bounce-back from barriers
    for (let y=0; y<ydim; y++) {
      for (let x=0; x<xdim; x++) {
        const i = index(x, y);
        const iq = i*q;
        if (barriers[i]) {
          for (let j=1;j<q;j++) {
            streamed[cIndex(x, y, 1, j)] = collided[iq + opp[j]];
          }
        }
      }
    }
  }

  public step(viscosity: number, barriers: boolean[]) {
    this.collide(viscosity);
    this.stream(barriers);
  }

  // Set all densities in a cell to their equilibrium values for a given velocity and density:
  public setEquilibrium(x: number, y: number, ux: number, uy: number, rho: number) {
    const { xdim, streamed } = this;
    const i = x + y*xdim;
    this._rho[i] = rho;

    const iq = i*q;
    const u2 =  1 - 1.5 * (ux * ux + uy * uy);
    for (let j = 0; j < q; j++) {
      const dir = cxs[j]*ux + cys[j]*uy;
      streamed[iq+j] =  weights[j] * rho * (u2 + 3 * dir + 4.5 * dir * dir);
    }
  }
}