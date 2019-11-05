export type Lattice = {
  dims: number[];
  vectors: number[][];
  weights: number[];
  opposites: number[];
};

function createVectors(size: number, weights: number[], q: number) {
  const wordsize = q*size;
  const buffer = new Float32Array(wordsize);
  for (let i = 0; i < wordsize; i+=q) {
    buffer.set(weights, i);
  }
  return buffer;
}

function collide(
  i: number, iq: number, viscosity: number,
  d: number, q: number,
  vectors: number[][], weights: number[],
  rho: Float32Array, streamed: Float32Array, collided: Float32Array
) {
  /* Calculate macroscopic quantities */

  let newrho = streamed[iq];       // macroscopic density
  const u = new Array(d).fill(0);  // macroscopic velocity components
  for (let j = 1; j < q; j++) {
    const uk = streamed[iq+j];
    newrho += uk;
    const v = vectors[j];
    for (let k=0; k<d; k++) {
      u[k] += v[k]*uk;
    }
  }

  // hack for stability
  if (newrho <= 0) newrho = 0.01;

  rho[i] = newrho;
  for (let k=0; k<d; k++) {
    u[k] /= newrho;
  }

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

  const u2 =  1 - 1.5 * u.reduce((a, n) => a + n*n, 0);
  streamed[iq] = omega * weights[0] * newrho * u2 + invomega * streamed[iq];
  for (let j = 1; j < q; j++) {
    const dir = vectors[j].reduce((a, c, i) => a + c*u[i], 0);
    const eq = weights[j] * newrho * (u2 + 3 * dir + 4.5 * dir * dir);
    collided[iq+j] = omega * eq + invomega * streamed[iq+j];
  }
}

function stream(
  i: number, iq: number,
  barriers: boolean[],
  q: number, opp: number[],
  destination: (i: number, j: number) => number,
  collided: Float32Array, streamed: Float32Array,
) {
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
  public rho: Float32Array;       // macroscopic density; cached for rendering
  private size: number;           // total number of lattice points
  private dims: number[];
  private vectors: number[][];
  private weights: number[];
  private opp: number[];
  private d: number;
  private q: number;

  constructor(lattice: Lattice) {
    const dims = this.dims = lattice.dims;
    const weights = this.weights = lattice.weights;
    const q = this.q = weights.length;
    const size = this.size = dims.reduce((a, n) => a * n, 1);
    this.streamed = createVectors(size, weights, q);
    this.collided = createVectors(size, weights, q);
    this.rho = new Float32Array(size);
    this.vectors = lattice.vectors;
    this.opp = lattice.opposites;
    this.d = dims.length;
  }


  public step(viscosity: number, barriers: boolean[]) {
    const { dims, size, vectors, weights, opp, d, q, rho, collided, streamed } = this;

    for (let i=0,iq=0; i<size; i++,iq+=q) {
      collide(i, iq, viscosity, d, q, vectors, weights, rho, streamed, collided);
    }

    const destination = (i: number, j: number) => {
      const v = vectors[j];
      let c = v[d-1];
      for (let k=d-2; k>=0; k--){
        c = c*dims[k] + v[k];
      }
      return q*((i + c + size) % size) + j;
    }
  
    for (let i=0,iq=0; i<size; i++,iq+=q) {
      stream(i, iq, barriers, q, opp, destination, collided, streamed);
    }
  }

  // Set all densities in a cell to their equilibrium values for a given velocity and density:
  public setEquilibrium(coords: number[], u: number[], rho: number) {
    const { dims, vectors, weights, d, q, streamed } = this;
  
    let i = coords[d-1];
    for (let j=d-2; j>=0; j--){
      i = i*dims[j] + coords[j];
    }

    this.rho[i] = rho;

    const iq = i*q;
    const u2 =  1 - 1.5 * u.reduce((a, n) => a + n*n, 0);
    for (let j = 0; j < q; j++) {
      const dir = vectors[j].reduce((a, c, i) => a + c*u[i], 0);
      streamed[iq+j] = weights[j] * rho * (u2 + 3 * dir + 4.5 * dir * dir);
    }
  }
}