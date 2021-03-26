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

function update_macros(max: number, macros: Float32Array, stationary: Float32Array, streamed: Float32Array) {
  for (let i=0,im=0; i<max; i++,im+=4) {
    let rho = stationary[i];  // macroscopic density
    let ux = 0, uy = 0, uz = 0; // macroscopic velocity components
    for (let j = 0, ij = i; j < q; j++,ij+=max) {
      const v = streamed[ij];
      const cj = c[j];
      rho += v;
      ux += cj[0]*v;
      uy += cj[1]*v;
      uz += cj[2]*v;
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
  max: number, omega: number, invomega: number,
  macros: Float32Array, collided: Float32Array, streamed: Float32Array,
) {
  
  /*const imax = max * q;
  for (let i = 0; i < imax; i++) {
    const j = (i / max)|0;
    const im = (i % max) << 2;*/

  const max4 = max << 2;
  for (let j=0,i=0; j < q; j++) {
    const [cx, cy, cz] = c[j];
    const omega_w = omega * weights[j];
    for (let im=0; im < max4; im+=4,i++) {
      const ux = macros[im+1];
      const uy = macros[im+2];
      const uz = macros[im+3];

      const u2 =  1 - 1.5 * (ux * ux + uy * uy + uz * uz);
      const dir = cx*ux + cy*uy + cz*uz;
      collided[i] = omega_w * macros[im] * (u2 + 3 * dir + 4.5 * dir * dir) +
                    invomega * streamed[i];
    }
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
  barriers: boolean[], streamed: Float32Array, collided: Float32Array,
) {
  /*const imax = max * q;
  for (let i = 0; i < imax; i++) {
    const j = (i / max)|0;
    const idx = i % max;*/
  for (let j=0,i=0; j < q; j++) {
    const [cx, cy, cz] = c[j];
    const oppindex = opp[j] * max;
    const srcindex = j * max;
    const shift = ((max<<1) - (cx-(cy-cz*ydim)*xdim))%max;
    for (let offset=shift; offset < max; offset++,i++) {
      const base = barriers[i] ? oppindex : srcindex;
      streamed[i] = collided[base + offset];
    }
    for (let offset=0; offset < shift; offset++,i++) {
      const base = barriers[i] ? oppindex : srcindex;
      streamed[i] = collided[base + offset];
    }
  }
}

export default class LatticeBoltzmann { 
  private stationary: Float32Array;       
  private streamed: Float32Array; // microscopic densities along each lattice direction
  private collided: Float32Array;
  public macros: Float32Array;      // velocity and density

  constructor(public xdim: number, public ydim: number, public zdim: number) {
    const size = xdim * ydim * zdim;

    const salloc = new ArrayBuffer(size * q << 2);
    for (let i = 0; i < q; i++)
      new Float32Array(salloc, i * size << 2, size).fill(weights[i]);

    const calloc = new ArrayBuffer(size * q << 2);
    for (let i = 0; i < q; i++)
      new Float32Array(calloc, i * size << 2, size);

    this.streamed = new Float32Array(salloc);
    this.collided = new Float32Array(calloc);
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
    collide(max, omega, invomega, macros, collided, streamed);
    stream(max, xdim, ydim, barriers, streamed, collided);
  }

  // Set all densities in a cell to their equilibrium values for a given velocity and density:
  public setEquilibrium(x: number, y: number, z: number, ux: number, uy: number, uz: number, rho: number) {
    const { xdim, ydim, zdim, streamed } = this;
    const max = xdim * ydim * zdim;
    const i = x+(y+z*ydim)*xdim;
    this.macros[i<<2] = rho;

    const u2 =  1 - 1.5 * (ux * ux + uy * uy + uz * uz);
    
    this.stationary[i] = w0 * rho * u2;
    for (let j = 0; j < q; j++) {
      const dir = c[j][0]*ux + c[j][1]*uy + c[j][2]*uz;
      streamed[j * max + i] = weights[j] * rho * (u2 + 3 * dir + 4.5 * dir * dir);
    }
  }
}