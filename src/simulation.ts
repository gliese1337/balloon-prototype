const w0 = 4/9;
const weights = [1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
const c: [number, number][] = [
  [ 0,  1 ], [  0, -1 ],
  [ 1,  0 ], [ -1,  0 ],
  [ 1,  1 ], [ -1, -1 ],
  [ 1, -1 ], [ -1,  1 ],
];

const opp = [ 1, 0, 3, 2, 5, 4, 7, 6 ];

const q = 8;

function update_macros(max: number, stationary: Float32Array, streamed: Float32Array[], macros: Float32Array) {
  for (let i = 0, im = 0; i < max; i++, im += 3) {
    let rho = stationary[i];  // macroscopic density
    let ux = 0, uy = 0; // macroscopic velocity components
    for (let j = 0; j < q; j++) {
      const v = streamed[j][i];
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
  [cx, cy]: [number, number],
  macros: Float32Array, streamed: Float32Array,
  collided: Float32Array,
) {
  for (let i = 0, im = 0; i < max; i++, im += 3) {
    const ux = macros[im+1];
    const uy = macros[im+2];

    const u2 =  1 - 1.5 * (ux * ux + uy * uy);
    const dir = cx*ux + cy*uy;
    collided[i] = omega_w * macros[im] * (u2 + 3 * dir + 4.5 * dir * dir) +
                          invomega * streamed[i];
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
  barriers: boolean[],
  collided: Float32Array,
  reflected: Float32Array,
  streamed: Float32Array,
) {
  for (let i = 0; i < max; i++) {
    const source = barriers[i] ? reflected : collided;
    streamed[i] = source[(i-(cx-cy*xdim)+max)%max];
  }
}

export default class LatticeBoltzmann {
  private stationary: Float32Array; // microscopic densities at 0 velocity
  private streamed: Float32Array[]; // microscopic densities along each lattice direction
  private collided: Float32Array[];
  public macros: Float32Array;    // macroscopic density & velocity

  constructor(public readonly xdim: number, public readonly ydim: number) {
    const size = xdim * ydim;

    this.stationary = new Float32Array(size).fill(w0);

    const plane_size = size * 4;
    const alloc_size = q * plane_size;
    const salloc = new ArrayBuffer(alloc_size);
    const calloc = new ArrayBuffer(alloc_size);
    const streamed: Float32Array[] = [];
    const collided: Float32Array[] = [];
    for (let i = 0, p = 0; i < q; i++, p += plane_size) {
      streamed.push(new Float32Array(salloc, p, size).fill(weights[i]));
      collided.push(new Float32Array(calloc, p, size));
    }

    this.streamed = streamed;
    this.collided = collided;
    this.macros = new Float32Array(size * 3);
  }

  public step(viscosity: number, barriers: boolean[]) {
    const { xdim, ydim, macros, collided, stationary, streamed } = this;
    const max = xdim * ydim;

    const tau = 3*viscosity + 0.5; // relaxation timescale
    const omega = 1 / tau;
    const invomega = 1 - omega;

    update_macros(max, stationary, streamed, macros);
    for (let j = 0; j < q; j++) {
      collide(max, omega * weights[j], invomega, c[j], macros, streamed[j], collided[j]);
    }

    update_static(max, omega * w0, invomega, macros, stationary);
    for (let j = 0; j < q; j++) {
      stream(xdim, max, c[j], barriers, collided[j], collided[opp[j]], streamed[j]);
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
    for (let j = 0, plane = 0; j < q; j++, plane += max) {
      const dir = c[j][0]*ux + c[j][1]*uy;
      streamed[j][i] = weights[j] * rho * (u2 + 3 * dir + 4.5 * dir * dir);
    }
  }
}