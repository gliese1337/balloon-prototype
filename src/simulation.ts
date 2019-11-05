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

function collide(i: number, iq: number, viscosity: number, rho: Float32Array, streamed: Float32Array, collided: Float32Array) {
  /* Calculate macroscopic quantities */

  let newrho = streamed[iq]; // macroscopic density
  let ux = 0, uy = 0;        // macroscopic velocity components
  for (let j = 1; j < q; j++) {
    const v = streamed[iq+j];
    newrho += v;
    ux += cxs[j]*v;
    uy += cys[j]*v;
  }

  // hack for stability
  if (newrho <= 0) newrho = 0.01;

  rho[i] = newrho;
  ux/=newrho;
  uy/=newrho;

  /* Calculate collisions */

  /*
  f_j^eq = w_j rho (1 + (e_j|u) / c_s^2 + (e_j|u)^2/(2 c_s^4) - u^2/(2 c_s^2))
  f`_j = omega f_j^eq + (1 - omega) f_j 

  c_s, the lattice speed of sound, is 1/sqrt(3)
  Thus, 1/c_s^2 = 3
        1/(2 c_s^4) = 4.5
        1/(2 c_s^2) = 1.5
  */

  const tau = 3*viscosity + 0.5; // relaxation timescale
  const omega = 1 / tau;
  const invomega = 1 - omega;

  const u2 =  1 - 1.5 * (ux * ux + uy * uy);
  streamed[iq] = omega * weights[0] * newrho * u2 + invomega * streamed[iq];
  for (let j = 1; j < q; j++) {
    const dir = cxs[j]*ux + cys[j]*uy;
    const eq = weights[j] * newrho * (u2 + 3 * dir + 4.5 * dir * dir);
    collided[iq+j] = omega * eq + invomega * streamed[iq+j];
  }
}

function stream(i: number, iq: number, barriers: boolean[], destination: (i: number, j: number) => number, collided: Float32Array, streamed: Float32Array) {
  if (barriers[i]) {
    // Handle bounce-back from barriers
    for (let j=1;j<q;j++) {
      streamed[destination(i, j)] = collided[iq + opp[j]];
    }
  } else {
    // Move particles along their velocity vector
    for (let j=1;j<q;j++) {
      streamed[destination(i, j)] = collided[iq + j];
    }
  }
}

export default class LatticeBoltzmann {        
  private streamed: Float32Array; // microscopic densities along each lattice direction
  private collided: Float32Array;
  public rho: Float32Array;    // macroscopic density; cached for rendering

  constructor(public readonly xdim: number, public readonly ydim: number) {
    const size = xdim * ydim;
    this.streamed = createVectors(size);
    this.collided = createVectors(size);
    this.rho = new Float32Array(size);
  }

  public step(viscosity: number, barriers: boolean[]) {
    const { xdim, ydim, rho, collided, streamed } = this;
    const max = xdim * ydim;

    for (let i=0,iq=0; i<max; i++,iq+=q) {
      collide(i, iq, viscosity, rho, streamed, collided);
    }

    const destination = (i: number, j: number) =>
                      q*((i+(cxs[j]+cys[j]*xdim)+max)%max)+j;

    for (let i=0,iq=0; i<max; i++,iq+=q) {
      stream(i, iq, barriers, destination, collided, streamed);
    }
  }

  // Set all densities in a cell to their equilibrium values for a given velocity and density:
  public setEquilibrium(x: number, y: number, ux: number, uy: number, rho: number) {
    const { xdim, streamed } = this;
    const i = x + y*xdim;
    this.rho[i] = rho;

    const iq = i*q;
    const u2 =  1 - 1.5 * (ux * ux + uy * uy);
    for (let j = 0; j < q; j++) {
      const dir = cxs[j]*ux + cys[j]*uy;
      streamed[iq+j] = weights[j] * rho * (u2 + 3 * dir + 4.5 * dir * dir);
    }
  }
}