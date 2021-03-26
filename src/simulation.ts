const weights = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
const c: [number, number][] = [
  [ 0,  0 ],
  [ 0,  1 ], [  0, -1 ],
  [ 1,  0 ], [ -1,  0 ],
  [ 1,  1 ], [ -1, -1 ],
  [ 1, -1 ], [ -1,  1 ],
]
const opp = [ 0, 2, 1, 4, 3, 6, 5, 8, 7];

const q = 9;

function createVectors(size: number) {
  const wordsize = q*size;
  const buffer = new Float32Array(wordsize);
  for (let i = 0; i < wordsize; i+=q) {
    buffer.set(weights, i);
  }
  return buffer;
}

function update_rho(max: number, streamed: Float32Array, rho: Float32Array) {
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

function update_velocities(max: number, rho: Float32Array, streamed: Float32Array, velocities: Float32Array) {
  for (let i=0,iv=0,iq=0; i<max; i++,iq+=q,iv+=2) {
    let ux = 0, uy = 0; // macroscopic velocity components
    for (let j = 1; j < q; j++) {
      const v = streamed[iq+j];
      ux += c[j][0]*v;
      uy += c[j][1]*v;
    }

    const newrho = rho[i];
    velocities[iv+0] = ux / newrho;
    velocities[iv+1] = uy / newrho;
  }
}

function collide(max: number, omega: number, invomega: number, rho: Float32Array, velocities: Float32Array, streamed: Float32Array, collided: Float32Array) {
  for (let j = 1; j < q; j++) {
    const [cx, cy] = c[j];
    const omega_w = omega * weights[j];
    for (let i=0,iv=0,iq=0; i<max; i++,iq+=q,iv+=2) {
      const ux = velocities[iv];
      const uy = velocities[iv+1];

      const u2 =  1 - 1.5 * (ux * ux + uy * uy);
      const dir = cx*ux + cy*uy;
      collided[iq+j] = omega_w * rho[i] * (u2 + 3 * dir + 4.5 * dir * dir) +
                       invomega * streamed[iq+j];
    }
  }
}

function stationary(max: number, omega: number, invomega: number, rho: Float32Array, velocities: Float32Array, streamed: Float32Array) {
  for (let j=1;j<q;j++) {
    for (let i=0,iv=0,iq=0; i<max; i++,iq+=q,iv+=2) {
      const ux = velocities[iv];
      const uy = velocities[iv+1];
      const u2 =  1 - 1.5 * (ux * ux + uy * uy);
      streamed[iq] = omega * weights[0] * rho[i] * u2 + invomega * streamed[iq];
    }
  }
}

function stream(xdim: number, max: number, barriers: boolean[], collided: Float32Array, streamed: Float32Array) {
  for (let j=1;j<q;j++) {
    const [cx, cy] = c[j];
    for (let i=0,iq=0; i<max; i++,iq+=q) {    
      if (barriers[i]) {
        streamed[iq + j] = collided[q*((i-(cx-cy*xdim)+max)%max) + opp[j]];
      } else {
        streamed[iq + j] = collided[q*((i-(cx-cy*xdim)+max)%max) + j];
      }
    }
  }
}

export default class LatticeBoltzmann {        
  private streamed: Float32Array; // microscopic densities along each lattice direction
  private collided: Float32Array;
  public rho: Float32Array;    // macroscopic density; cached for rendering
  public velocities: Float32Array; // macroscopic velocity

  constructor(public readonly xdim: number, public readonly ydim: number) {
    const size = xdim * ydim;
    this.streamed = createVectors(size);
    this.collided = createVectors(size);
    this.rho = new Float32Array(size);
    this.velocities = new Float32Array(size*2);
  }

  public step(viscosity: number, barriers: boolean[]) {
    const { xdim, ydim, rho, velocities, collided, streamed } = this;
    const max = xdim * ydim;

    const tau = 3*viscosity + 0.5; // relaxation timescale
    const omega = 1 / tau;
    const invomega = 1 - omega;

    update_rho(max, streamed, rho);
    update_velocities(max, rho, streamed, velocities);
    collide(max, omega, invomega, rho, velocities, streamed, collided);
    stationary(max, omega, invomega, rho, velocities, streamed);
    stream(xdim, max, barriers, collided, streamed);
  }

  // Set all densities in a cell to their equilibrium values for a given velocity and density:
  public setEquilibrium(x: number, y: number, ux: number, uy: number, rho: number) {
    const { xdim, streamed } = this;
    const i = x + y*xdim;
    this.rho[i] = rho;

    const iq = i*q;
    const u2 =  1 - 1.5 * (ux * ux + uy * uy);
    for (let j = 0; j < q; j++) {
      const dir = c[j][0]*ux + c[j][1]*uy;
      streamed[iq+j] = weights[j] * rho * (u2 + 3 * dir + 4.5 * dir * dir);
    }
  }
}