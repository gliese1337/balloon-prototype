const w0 = 4/9;
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

function update_macros(max: number, stationary: Float32Array, streamed: Float32Array, macros: Float32Array) {
  for (let i = 0, im = 0; i < max; i++, im += 3) {
    let rho = stationary[i];  // macroscopic density
    let ux = 0, uy = 0; // macroscopic velocity components
    for (let j = 1, plane = max; j < q; j++, plane += max) {
      const v = streamed[plane + i];
      rho += v;
      ux += c[j][0]*v;
      uy += c[j][1]*v;
    }

    // hack for stability
    rho = rho <= 0 ? 0.01 : rho;

    macros[im] = rho;
    macros[im + 1] = ux / rho;
    macros[im + 2] = uy / rho;
  }
}

function collide(
  max: number, omega_w: number, invomega: number,
  [cx, cy]: [number, number], plane: number,
  macros: Float32Array, streamed: Float32Array,
  collided: Float32Array,
) {
  for (let i = 0, im = 0; i < max; i++, im += 3) {
    const ux = macros[im+1];
    const uy = macros[im+2];

    const u2 =  1 - 1.5 * (ux * ux + uy * uy);
    const dir = cx*ux + cy*uy;
    collided[plane + i] = omega_w * macros[im] * (u2 + 3 * dir + 4.5 * dir * dir) +
                          invomega * streamed[plane + i];
  }
}

function update_static(max: number, omega_w: number, invomega: number, macros: Float32Array, stationary: Float32Array) {
  for (let i = 0, im = 0; i < max; i++, im += 3) {
    const ux = macros[im + 1];
    const uy = macros[im + 2];
    const u2 =  1 - 1.5 * (ux * ux + uy * uy);
    stationary[i] = omega_w * macros[im] * u2 + invomega * stationary[i];
  }
}

function stream(
  xdim: number, max: number,
  [cx, cy]: [number, number],
  plane: number, opp_plane: number,
  barriers: boolean[], collided: Float32Array,
  streamed: Float32Array,
) {
  for (let i = 0; i < max; i++) {
    const source = barriers[i] ? opp_plane : plane;
    streamed[plane + i] = collided[source + ((i-(cx-cy*xdim)+max)%max)];
  }
}

export default class LatticeBoltzmann {
  private stationary: Float32Array;
  private streamed: Float32Array; // microscopic densities along each lattice direction
  private collided: Float32Array;
  public macros: Float32Array;    // macroscopic density & velocity

  constructor(public readonly xdim: number, public readonly ydim: number) {
    const size = xdim * ydim;

    this.stationary = new Float32Array(size).fill(w0);

    const salloc = new ArrayBuffer(size * q << 2);
    for (let i = 0; i < q; i++)
      new Float32Array(salloc, i * size << 2, size).fill(weights[i]);

    this.streamed = new Float32Array(salloc);
    this.collided = new Float32Array(size * q);
    this.macros = new Float32Array(size * 3);
  }

  public step(viscosity: number, barriers: boolean[]) {
    const { xdim, ydim, macros, collided, stationary, streamed } = this;
    const max = xdim * ydim;

    const tau = 3*viscosity + 0.5; // relaxation timescale
    const omega = 1 / tau;
    const invomega = 1 - omega;

    update_macros(max, stationary, streamed, macros);
    for (let j = 1, plane = max; j < q; j++, plane += max) {
      const omega_w = omega * weights[j];
      collide(max, omega_w, invomega, c[j], plane, macros, streamed, collided);
    }

    update_static(max, omega * w0, invomega, macros, stationary);
    for (let j = 1, plane = max; j < q; j++, plane += max) {
      const opp_plane = max * opp[j];
      stream(xdim, max, c[j], plane, opp_plane, barriers, collided, streamed);
    }
  }

  // Set all densities in a cell to their equilibrium values for a given velocity and density:
  public setEquilibrium(x: number, y: number, ux: number, uy: number, rho: number) {
    const { xdim, ydim, streamed } = this;
    const max = xdim * ydim;
    const i = x + y*xdim;
    this.macros[i * 3] = rho;

    const u2 =  1 - 1.5 * (ux * ux + uy * uy);
    
    this.stationary[i] = w0 * rho * u2;
    for (let j = 1, plane = max; j < q; j++, plane += max) {
      const dir = c[j][0]*ux + c[j][1]*uy;
      streamed[plane + i] = weights[j] * rho * (u2 + 3 * dir + 4.5 * dir * dir);
    }
  }
}