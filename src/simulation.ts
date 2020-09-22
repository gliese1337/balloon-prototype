const weights = [1/3, 1/18, 1/18, 1/18, 1/18, 1/18, 1/18, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36];
const cxs =     [  0,    1,   -1,    0,    0,    0,    0,    1,   -1,   -1,    1,    1,   -1,    1,   -1,    0,    0,    0,    0];
const cys =     [  0,    0,    0,    1,   -1,    0,    0,    1,   -1,    1,   -1,    0,    0,    0,    0,    1,   -1,   -1,    1];
const czs =     [  0,    0,    0,    0,    0,    1,   -1,    0,    0,    0,    0,    1,   -1,   -1,    1,    1,   -1,    1,   -1];
const opp =     [  0,    2,    1,    4,    3,    6,    5,    8,    7,   10,    9,   12,   11,   14,   13,   16,   15,   18,   17];
const q = 19;

function createVectors(size: number) {
  const wordsize = q*size;
  const buffer = new Float32Array(wordsize);
  for (let i = 0; i < wordsize; i+=q) {
    buffer.set(weights, i);
  }
  return buffer;
}

function update_macros(max: number, macros: Float32Array, streamed: Float32Array) {
  for (let i=0,im=0,iq=0; i<max; i++,im+=4,iq+=q) {
    let newrho = streamed[iq];  // macroscopic density
    let ux = 0, uy = 0, uz = 0; // macroscopic velocity components
    for (let j = 1; j < q; j++) {
      const v = streamed[iq+j];
      newrho += v;
      ux += cxs[j]*v;
      uy += cys[j]*v;
      uz += czs[j]*v;
    }

    // hack for stability
    newrho = newrho <= 0 ? 0.01 : newrho;

    macros[im] = newrho;
    macros[im+1] = ux / newrho;
    macros[im+2] = uy / newrho;
    macros[im+3] = uz / newrho;
  }
}

function collide(max: number, viscosity: number, macros: Float32Array, collided: Float32Array, streamed: Float32Array) {
  const tau = 3*viscosity + 0.5; // relaxation timescale
  const omega = 1 / tau;
  const invomega = 1 - omega;

  for (let i=0,im=0,iq=0; i<max; i++,im+=4,iq+=q) {
    const newrho = macros[im];
    const ux = macros[im+1];
    const uy = macros[im+2];
    const uz = macros[im+3];

    /* Calculate collisions */

    /*
    f_j^eq = w_j rho (1 + (e_j|u) / c_s^2 + (e_j|u)^2/(2 c_s^4) - u^2/(2 c_s^2))
    f`_j = omega f_j^eq + (1 - omega) f_j 

    c_s, the lattice speed of sound, is 1/sqrt(3)
    Thus, 1/c_s^2 = 3
          1/(2 c_s^4) = 4.5
          1/(2 c_s^2) = 1.5
    */

    const u2 =  1 - 1.5 * (ux * ux + uy * uy + uz * uz);
    for (let j = 1; j < q; j++) {
      const dir = cxs[j]*ux + cys[j]*uy + czs[j]*uz;
      const eq = weights[j] * newrho * (u2 + 3 * dir + 4.5 * dir * dir);
      collided[iq+j] = omega * eq + invomega * streamed[iq+j];
    }
  }
}

function stream(xdim: number, ydim: number, max: number, viscosity: number, barriers: boolean[], macros: Float32Array, collided: Float32Array, streamed: Float32Array) {
  const tau = 3*viscosity + 0.5; // relaxation timescale
  const omega = 1 / tau;
  const invomega = 1 - omega;

  for (let i=0,im=0,iq=0; i<max; i++,im+=4,iq+=q) {
    const ux = macros[im+1];
    const uy = macros[im+2];
    const uz = macros[im+3];

    const u2 =  1 - 1.5 * (ux * ux + uy * uy + uz * uz);
    streamed[iq] = omega * weights[0] * macros[im] * u2 + invomega * streamed[iq];
    
    if (barriers[i]) {
      // Handle bounce-back from barriers
      for (let j=1;j<q;j++) {
        streamed[q*((i+(cxs[j]+(cys[j]+czs[j]*ydim)*xdim)+max)%max)+j] = collided[iq + opp[j]];
      }
    } else {
      // Move particles along their velocity vector
      for (let j=1;j<q;j++) {
        streamed[q*((i+(cxs[j]+(cys[j]+czs[j]*ydim)*xdim)+max)%max)+j] = collided[iq + j];
      }
    }
  }
}

export default class LatticeBoltzmann {        
  private streamed: Float32Array; // microscopic densities along each lattice direction
  private collided: Float32Array;
  public macros: Float32Array;        // velocity and density

  constructor(public xdim: number, public ydim: number, public zdim: number) {
    const size = xdim * ydim * zdim;
    this.streamed = createVectors(size);
    this.collided = createVectors(size);
    this.macros = new Float32Array(size*4);
  }

  public step(viscosity: number, barriers: boolean[]) {
    const { xdim, ydim, zdim, macros, collided, streamed } = this;
    const max = xdim * ydim * zdim;
    update_macros(max, macros, streamed);
    collide(max, viscosity, macros, collided, streamed);
    stream(xdim, ydim, max, viscosity, barriers, macros, collided, streamed);
  }

  // Set all densities in a cell to their equilibrium values for a given velocity and density:
  public setEquilibrium(x: number, y: number, z: number, ux: number, uy: number, uz: number, rho: number) {
    const { xdim, ydim, streamed } = this;
    const i = x+(y+z*ydim)*xdim;
    this.macros[i<<2] = rho;

    const iq = i*q;
    const u2 =  1 - 1.5 * (ux * ux + uy * uy + uz * uz);
    for (let j = 0; j < q; j++) {
      const dir = cxs[j]*ux + cys[j]*uy;
      streamed[iq+j] = weights[j] * rho * (u2 + 3 * dir + 4.5 * dir * dir);
    }
  }
}