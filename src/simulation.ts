const weights = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
const c: [number, number][] = [
  [ 0,  0 ],
  [ 0,  1 ], [  0, -1 ],
  [ 1,  0 ], [ -1,  0 ],
  [ 1,  1 ], [ -1, -1 ],
  [ 1, -1 ], [ -1,  1 ],
];

const opp = [ 0, 2, 1, 4, 3, 6, 5, 8, 7 ];

const q = 9;

function update_rho(max: number, streamed: Float32Array, rho: Float32Array) {
  for (let i = 0; i < max; i++) {
    let newrho = streamed[i];  // macroscopic density
    for (let j = 1, plane = max; j < q; j++, plane += max) {
      const v = streamed[plane + i];
      newrho += v;
    }

    // hack for stability
    rho[i] = newrho <= 0 ? 0.01 : newrho;
  }
}

function update_velocities(max: number, rho: Float32Array, streamed: Float32Array, velocities: Float32Array) {
  for (let i = 0,iv = 0; i < max; i++, iv += 2) {
    let ux = 0, uy = 0; // macroscopic velocity components
    for (let j = 1, plane = max; j < q; j++, plane += max) {
      const v = streamed[plane + i];
      ux += c[j][0]*v;
      uy += c[j][1]*v;
    }

    const newrho = rho[i];
    velocities[iv+0] = ux / newrho;
    velocities[iv+1] = uy / newrho;
  }
}

function collide(max: number, omega: number, invomega: number, rho: Float32Array, velocities: Float32Array, streamed: Float32Array, collided: Float32Array) {
  for (let j = 1, plane = max; j < q; j++, plane += max) {
    const [cx, cy] = c[j];
    const omega_w = omega * weights[j];
    for (let i = 0, iv = 0; i < max; i++, iv += 2) {
      const ux = velocities[iv];
      const uy = velocities[iv+1];

      const u2 =  1 - 1.5 * (ux * ux + uy * uy);
      const dir = cx*ux + cy*uy;
      collided[plane + i] = omega_w * rho[i] * (u2 + 3 * dir + 4.5 * dir * dir) +
                       invomega * streamed[plane + i];
    }
  }
}

function stationary(max: number, omega: number, invomega: number, rho: Float32Array, velocities: Float32Array, streamed: Float32Array) {
  for (let i = 0, iv = 0; i < max; i++, iv += 2) {
    const ux = velocities[iv];
    const uy = velocities[iv+1];
    const u2 =  1 - 1.5 * (ux * ux + uy * uy);
    streamed[i] = omega * weights[0] * rho[i] * u2 + invomega * streamed[i];
  }
}

function stream(xdim: number, max: number, barriers: boolean[], collided: Float32Array, streamed: Float32Array) {
  for (let j = 1, plane = max; j < q; j++, plane += max) {
    const [cx, cy] = c[j];
    const opp_plane = max * opp[j];
    for (let i=0; i<max; i++) {
      const source = barriers[i] ? opp_plane : plane;
      streamed[plane + i] = collided[source + ((i-(cx-cy*xdim)+max)%max)];
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

    const salloc = new ArrayBuffer(size * q << 2);
    for (let i = 0; i < q; i++)
      new Float32Array(salloc, i * size << 2, size).fill(weights[i]);

    this.streamed = new Float32Array(salloc);
    this.collided = new Float32Array(size * q);
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
    const { xdim, ydim, streamed } = this;
    const max = xdim * ydim;
    const i = x + y*xdim;
    this.rho[i] = rho;

    const u2 =  1 - 1.5 * (ux * ux + uy * uy);
    for (let j = 0, plane = 0; j < q; j++, plane += max) {
      const dir = c[j][0]*ux + c[j][1]*uy;
      streamed[plane + i] = weights[j] * rho * (u2 + 3 * dir + 4.5 * dir * dir);
    }
  }
}