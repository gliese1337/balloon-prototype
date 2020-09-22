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

function update_rho(max: number, rho: Float32Array, streamed: Float32Array) {
  for (let i=0,iq=0; i<max; i++,iq+=q) {
    let newrho = streamed[iq];  // macroscopic density
    for (let j = 1; j < q; j++) {
      const v = streamed[iq+j];
      newrho += v;
    }

    // hack for stability
    rho[i] = newrho <= 0 ? 0.01 : newrho;
  }
}

function update_velocities(max: number, rho: Float32Array, velocities: Float32Array, streamed: Float32Array) {
  for (let i=0,iv=0,iq=0; i<max; i++,iq+=q,iv+=3) {
    let ux = 0, uy = 0, uz = 0; // macroscopic velocity components
    for (let j = 1; j < q; j++) {
      const v = streamed[iq+j];
      ux += cxs[j]*v;
      uy += cys[j]*v;
      uz += czs[j]*v;
    }

    const newrho = rho[i];
    velocities[iv+0] = ux / newrho;
    velocities[iv+1] = uy / newrho;
    velocities[iv+2] = uz / newrho;
  }
}

function collide(max: number, viscosity: number, rho: Float32Array, velocities: Float32Array, collided: Float32Array, streamed: Float32Array) {
  const tau = 3*viscosity + 0.5; // relaxation timescale
  const omega = 1 / tau;
  const invomega = 1 - omega;

  for (let i=0,iv=0,iq=0; i<max; i++,iq+=q,iv+=3) {
    const ux = velocities[iv];
    const uy = velocities[iv+1];
    const uz = velocities[iv+2];
    const newrho = rho[i];

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

function stream(xdim: number, ydim: number, max: number, viscosity: number, rho: Float32Array, barriers: boolean[], velocities: Float32Array, collided: Float32Array, streamed: Float32Array) {
  const tau = 3*viscosity + 0.5; // relaxation timescale
  const omega = 1 / tau;
  const invomega = 1 - omega;

  for (let i=0,iv=0,iq=0; i<max; i++,iq+=q,iv+=3) {
    const ux = velocities[iv];
    const uy = velocities[iv+1];
    const uz = velocities[iv+2];

    const u2 =  1 - 1.5 * (ux * ux + uy * uy + uz * uz);
    streamed[iq] = omega * weights[0] * rho[i] * u2 + invomega * streamed[iq];
    
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
  public rho: Float32Array;        // macroscopic density; cached for rendering
  public velocities: Float32Array; // macroscopic velocity

  constructor(public xdim: number, public ydim: number, public zdim: number) {
    const size = xdim * ydim * zdim;
    this.streamed = createVectors(size);
    this.collided = createVectors(size);
    this.rho = new Float32Array(size);
    this.velocities = new Float32Array(size*3);
  }

  public step(viscosity: number, barriers: boolean[]) {
    const { xdim, ydim, zdim, rho, velocities, collided, streamed } = this;
    const max = xdim * ydim * zdim;
    update_rho(max, rho, streamed);
    update_velocities(max, rho, velocities, streamed);
    collide(max, viscosity, rho, velocities, collided, streamed);
    stream(xdim, ydim, max, viscosity, rho, barriers, velocities, collided, streamed);
  }

  // Set all densities in a cell to their equilibrium values for a given velocity and density:
  public setEquilibrium(x: number, y: number, z: number, ux: number, uy: number, uz: number, rho: number) {
    const { xdim, ydim, streamed } = this;
    const i = x+(y+z*ydim)*xdim;
    this.rho[i] = rho;

    const iq = i*q;
    const u2 =  1 - 1.5 * (ux * ux + uy * uy + uz * uz);
    for (let j = 0; j < q; j++) {
      const dir = cxs[j]*ux + cys[j]*uy;
      streamed[iq+j] = weights[j] * rho * (u2 + 3 * dir + 4.5 * dir * dir);
    }
  }
}