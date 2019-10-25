const four9ths = 4.0 / 9.0;
const one9th = 1.0 / 9.0;
const one36th = 1.0 / 36.0;

export default class LatticeBoltzmann {
  private n0: Float32Array;            // microscopic densities along each lattice direction
  private nN: Float32Array;
  private nS: Float32Array;
  private nE: Float32Array;
  private nW: Float32Array;
  private nNE: Float32Array;
  private nSE: Float32Array;
  private nNW: Float32Array;
  private nSW: Float32Array;
  
  public rho: Float32Array;           // macroscopic density

  constructor(public xdim: number, public ydim: number) {
    // Create the arrays of fluid particle densities, etc. (using 1D arrays for speed):
    // To index into these arrays, use x + y*xdim, traversing rows first and then columns.
    const size = 4*xdim * ydim;

    const buffer = new ArrayBuffer(40*size);
    const n0 =  this.n0 =  new Float32Array(buffer, 0, size);
    const nN =  this.nN =  new Float32Array(buffer, 4*size, size);
    const nS =  this.nS =  new Float32Array(buffer, 8*size, size);
    const nE =  this.nE =  new Float32Array(buffer, 12*size, size);
    const nW =  this.nW =  new Float32Array(buffer, 16*size, size);
    const nNE = this.nNE = new Float32Array(buffer, 20*size, size);
    const nSE = this.nSE = new Float32Array(buffer, 24*size, size);
    const nNW = this.nNW = new Float32Array(buffer, 28*size, size);
    const nSW = this.nSW = new Float32Array(buffer, 32*size, size);
    const rho = this.rho = new Float32Array(buffer, 36*size, size);

    for (let i=0; i<size; i++) {
      n0[i]  = four9ths;
      nE[i]  =   one9th;
      nW[i]  =   one9th;
      nN[i]  =   one9th;
      nS[i]  =   one9th;
      nNE[i] =  one36th;
      nSE[i] =  one36th;
      nNW[i] =  one36th;
      nSW[i] =  one36th;
      rho[i] = 1;
    }
  }

  // Collide particles within each cell (here's the physics!):
  public collide(viscosity: number) {
    const { xdim, ydim, rho, n0, nN, nNW, nE, nNE, nS, nSE, nW, nSW } = this;
    const omega = 1 / (3*viscosity + 0.5);        // reciprocal of relaxation time
    const oneminuso = 1 - omega;
    const max = xdim * ydim;
    for (let i=0; i < max; i++) {
      const n0i = n0[i];
      const nNi = nN[i];
      const nSi = nS[i];
      const nEi = nE[i];
      const nWi = nW[i];
      const nNEi = nNE[i];
      const nNWi = nNW[i];
      const nSEi = nSE[i];
      const nSWi = nSW[i];

      const newrho = n0i + nNi + nSi + nEi + nWi + nNWi + nNEi + nSWi + nSEi;
      rho[i] = newrho;

      // velocity components
      const ux = (nEi + nNEi + nSEi - nWi - nNWi - nSWi) / newrho;
      const uy = (nNi + nNEi + nNWi - nSi - nSEi - nSWi) / newrho;

      // common subexpressions
      const ux3 = 3 * ux;
      const uy3 = 3 * uy;
      const ux2 = ux * ux;
      const ux245 = 4.5 * ux2;
      const uy2 = uy * uy;
      const uy245 = 4.5 * uy2;
      const uxuy2 = 2 * ux * uy;
      const uxuy245 = 4.5 * uxuy2;
      const u2 = ux2 + uy2;
      const u245 = ux245 + uy245;
      const u215 = 1 - 1.5 * u2;
      const one9thrho = one9th * newrho * omega;
      const one36thrho = one36th * newrho * omega;
      n0[i]  = omega * four9ths*newrho * u215 + n0i*oneminuso;
      nE[i]  =  one9thrho * (u215 + ux3       + ux245       ) + nEi*oneminuso;
      nW[i]  =  one9thrho * (u215 - ux3       + ux245       ) + nWi*oneminuso;
      nN[i]  =  one9thrho * (u215 + uy3       + uy245       ) + nNi*oneminuso;
      nS[i]  =  one9thrho * (u215 - uy3       + uy245       ) + nSi*oneminuso;
      nNE[i] = one36thrho * (u215 + ux3 + uy3 + u245+uxuy245) + nNEi*oneminuso;
      nSE[i] = one36thrho * (u215 + ux3 - uy3 + u245-uxuy245) + nSEi*oneminuso;
      nNW[i] = one36thrho * (u215 - ux3 + uy3 + u245-uxuy245) + nNWi*oneminuso;
      nSW[i] = one36thrho * (u215 - ux3 - uy3 + u245+uxuy245) + nSWi*oneminuso;
    }

    for (let y=xdim-1; y<max; y+=xdim) {
      // at right end, copy left-flowing densities from next column to the left
      nW[y] = nW[y-1];
      nNW[y] = nNW[y-1];
      nSW[y] = nSW[y-1];
    }
  }

  // Move particles along their directions of motion:
  public stream(barriers: boolean[]) {
    const { xdim, ydim, nN, nNW, nE, nNE, nS, nSE, nW, nSW } = this;

    for (let y=(ydim-2)*xdim; y>0; y-=xdim) {
      const xlimit = xdim + y - 1;
      for (let x=1+y; x<xlimit; x++) {            // first start in NW corner...
        nN[x] = nN[x-xdim];            // move the north-moving particles
        nNW[x] = nNW[x+1-xdim];        // and the northwest-moving particles
      }
      for (let x=xdim-2+y; x>y; x--) {            // now start in NE corner...
        nE[x] = nE[x-1];            // move the east-moving particles
        nNE[x] = nNE[x-1-xdim];        // and the northeast-moving particles
      }
    }
    const ylimit = (ydim-1)*xdim;
    for (let y=xdim; y<ylimit; y+=xdim) {
      for (let x=xdim-2+y; x>y; x--) {             // now start in SE corner...
        nS[x] = nS[x+xdim];            // move the south-moving particles
        nSE[x] = nSE[x-1+xdim];        // and the southeast-moving particles
      }
      const xlimit = xdim + y - 1;
      for (let x=1+y; x<xlimit; x++) { // now start in the SW corner...
        nW[x] = nW[x+1];            // move the west-moving particles
        nSW[x] = nSW[x+1+xdim];        // and the southwest-moving particles
      }
    }
    for (let y=xdim; y<ylimit; y+=xdim) {                // Now handle bounce-back from barriers
      for (let x=1+y; x<y+xdim-1; x++) {
        if (barriers[x]) {
            nE[x+1] = nW[x];
            nW[x-1] = nE[x];
            nN[x+xdim] = nS[x];
            nS[x-xdim] = nN[x];
            nNE[x+1+xdim] = nSW[x];
            nNW[x-1+xdim] = nSE[x];
            nSE[x+1-xdim] = nNW[x];
            nSW[x-1-xdim] = nNE[x];
            // Force on the barrier:
            // barrierFx += nE[index] + nNE[index] + nSE[index] - nW[index] - nNW[index] - nSW[index];
            // barrierFy += nN[index] + nNE[index] + nNW[index] - nS[index] - nSE[index] - nSW[index];
        }
      }
    }
  }

  // Set all densities in a cell to their equilibrium values for a given velocity and density:
  public setEquilibrium(x: number, y: number, ux: number, uy: number, rho: number) {
    const { xdim, n0, nN, nNW, nE, nNE, nS, nSE, nW, nSW } = this;
    const i = x + y*xdim;
    this.rho[i] = rho;    const ux3 = 3 * ux;
    const uy3 = 3 * uy;
    const ux2 = ux * ux;
    const ux245 = 4.5 * ux2;
    const uy2 = uy * uy;
    const uy245 = 4.5 * uy2;
    const uxuy2 = 2 * ux * uy;
    const uxuy245 = 4.5 * uxuy2;
    const u2 = ux2 + uy2;
    const u245 = ux245 + uy245;
    const u215 = 1 - 1.5 * u2;
    const one9thrho = one9th * rho;
    const one36thrho = one36th * rho;
    n0[i]  = four9ths*rho * u215;
    nE[i]  =  one9thrho * (u215 + ux3       + ux245       );
    nW[i]  =  one9thrho * (u215 - ux3       + ux245       );
    nN[i]  =  one9thrho * (u215 + uy3       + uy245       );
    nS[i]  =  one9thrho * (u215 - uy3       + uy245       );
    nNE[i] = one36thrho * (u215 + ux3 + uy3 + u245+uxuy245);
    nSE[i] = one36thrho * (u215 + ux3 - uy3 + u245-uxuy245);
    nNW[i] = one36thrho * (u215 - ux3 + uy3 + u245-uxuy245);
    nSW[i] = one36thrho * (u215 - ux3 - uy3 + u245+uxuy245);
  }
}