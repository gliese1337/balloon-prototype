const w0 = 1/3;
const weights = [1/18, 1/18, 1/18, 1/18, 1/18, 1/18, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36];
const c: [number, number, number][] = [
  [ 1, 0, 0], [-1, 0, 0],
  [ 0, 1, 0], [0, -1, 0],
  [ 0, 0, 1], [ 0, 0, -1],
  [ 1, 1, 0], [-1,-1, 0], [-1, 1, 0], [ 1,-1, 0],
  [ 1, 0, 1], [-1, 0,-1], [ 1, 0,-1], [-1, 0, 1],
  [ 0, 1, 1], [ 0,-1,-1], [ 0,-1, 1], [ 0, 1,-1],
];
const opp = [1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 17, 16];
const q = 18;

function update_macros(max: number, macros: Float32Array, stationary: Float32Array, streamed: Float32Array[]) {
  for (let i=0,im=0; i<max; i++,im+=4) {
    let rho = stationary[i];  // macroscopic density
    let ux = 0, uy = 0, uz = 0; // macroscopic velocity components
    for (let j = 0; j < q; j++) {
      const v = streamed[j][i];
      rho += v;
      ux += c[j][0]*v;
      uy += c[j][1]*v;
      uz += c[j][2]*v;
    }

    // hack for stability
    rho = rho <= 0 ? 0.01 : rho;

    macros[im] = rho;
    macros[im+1] = ux / rho;
    macros[im+2] = uy / rho;
    macros[im+3] = uz / rho;
  }
}

function collide(
  max: number,
  omega_w: number, invomega: number,
  c: [number, number, number], macros: Float32Array,
  collided: Float32Array, streamed: Float32Array,
) {
  const [cx, cy, cz] = c;
  for (let i=0,im=0; i < max; i++,im+=4) {
    const ux = macros[im+1];
    const uy = macros[im+2];
    const uz = macros[im+3];

    const u2 =  1 - 1.5 * (ux * ux + uy * uy + uz * uz);
    const dir = cx*ux + cy*uy + cz*uz;
    collided[i] = omega_w * macros[im] * (u2 + 3 * dir + 4.5 * dir * dir) +
                  invomega * streamed[i];
  }
}

function update_static(max: number, omega_w: number, invomega: number, macros: Float32Array, stationary: Float32Array) {
  for (let i=0,im=0; i<max; i++,im+=4) {
    const ux = macros[im+1];
    const uy = macros[im+2];
    const uz = macros[im+3];

    const u2 =  1 - 1.5 * (ux * ux + uy * uy + uz * uz);
    stationary[i] = omega_w * macros[im] * u2 + invomega * stationary[i];
  }
}

function stream(
  max: number, xdim: number, ydim: number,
  c: [number, number, number],
  barriers: boolean[], streamed: Float32Array,
  collided: Float32Array, opposite: Float32Array,
) {
  const [cx, cy, cz] = c;
  const offset = cx-(cy-cz*ydim)*xdim;
  for (let i=0,im=0; i<max; i++,im+=4) {
    streamed[i] = (barriers[i] ? opposite : collided)[(i-offset+max)%max];
  }
}

export default class LatticeBoltzmann { 
  private stationary: Float32Array;       
  private streamed: Float32Array[]; // microscopic densities along each lattice direction
  private collided: Float32Array[];
  public macros: Float32Array;      // velocity and density

  constructor(public xdim: number, public ydim: number, public zdim: number) {
    const size = xdim * ydim * zdim;

    const salloc = new ArrayBuffer(size * q << 2);
    const streamed: Float32Array[] = [];
    for (let i = 0; i < q; i++)
      streamed.push(new Float32Array(salloc, i * size << 2, size).fill(weights[i]));

    const calloc = new ArrayBuffer(size * q << 2);
    const collided: Float32Array[] = [];
    for (let i = 0; i < q; i++)
      collided.push(new Float32Array(calloc, i * size << 2, size));

    this.streamed = streamed;
    this.collided = collided;
    this.macros = new Float32Array(size*4);
    this.stationary = new Float32Array(size).fill(w0);
  }

  public step(viscosity: number, barriers: boolean[]) {
    const { xdim, ydim, zdim, macros, collided, streamed, stationary } = this;
    const max = xdim * ydim * zdim;
    
    const tau = 3*viscosity + 0.5; // relaxation timescale
    const omega = 1 / tau;
    const invomega = 1 - omega;

    update_macros(max, macros, stationary, streamed);
    update_static(max, omega * w0, invomega, macros, stationary);
    for (let j = 0; j < q; j++)
      collide(max, omega * weights[j], invomega, c[j], macros, collided[j], streamed[j]);
    for (let j = 0; j < q; j++)
      stream(max, xdim, ydim, c[j], barriers, streamed[j], collided[j], collided[opp[j]]);
  }

  // Set all densities in a cell to their equilibrium values for a given velocity and density:
  public setEquilibrium(x: number, y: number, z: number, ux: number, uy: number, uz: number, rho: number) {
    const { xdim, ydim, streamed } = this;
    const i = x+(y+z*ydim)*xdim;
    this.macros[i<<2] = rho;

    const u2 =  1 - 1.5 * (ux * ux + uy * uy + uz * uz);
    
    this.stationary[i] = w0 * rho * u2;
    for (let j = 0; j < q; j++) {
      const dir = c[j][0]*ux + c[j][1]*uy + c[j][2]*uz;
      streamed[j][i] = weights[j] * rho * (u2 + 3 * dir + 4.5 * dir * dir);
    }
  }
}