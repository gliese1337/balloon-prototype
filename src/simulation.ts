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

function stream(
  i: number, iq: number,
  barriers: boolean[],
  q: number, opp: number[],
  destination: Function,
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

function dot(v: number[] | Float32Array, k: number[] | Float32Array, l: number) {
  let d = v[0]*k[0];
  for (let i=1;i<l;i++) d += v[i]*k[i];
  return d;
}

function destination(d: number, q: number, dims: number[], size: number, vectors: number[][]){
  let term = `v[${d-1}]*${dims[d-2]}`;
  for (let i=d-2;i>0;i--) {
    term = `(v[${i}]+${term})*${dims[i-1]}`
  }
  return new Function(
    'vectors', 'i', 'j',
    `const v = vectors[j]; return ${q}*((i+(v[0] + ${term})+${size})%${size})+j;`,
  ).bind(null, vectors);
}

function expandVectorLoop(d: number, q: number, vectors: number[][], weights: number[]) {
  const body = [];
  for (let j = 1; j < q; j++) {
    const v = vectors[j];
    const terms = [];
    for (let k=0;k<d;k++) {
      if (v[k] === 0) continue;
      terms.push(`${v[k] < 0 ? '-' : '+'}u${k}`);
    }
    
    const direxp = terms.join('');

    body.push(`dir = ${ direxp[0] === '+' ? direxp.substr(1) : direxp }; collided[iq+${j}] = omega * (${weights[j]} * newrho * (usqr + 3 * dir + 4.5 * dir * dir)) + invomega * streamed[iq+${j}]`);
  }

  return body.join('\n');
}

function collide(d: number, q: number, vectors: number[][], weights: number[]) {
  const body = `
    let newrho = streamed[iq];
    let ${ Array.from({ length: d }, (_, i) => `u${i} = 0`).join(', ') };
    for (let j = 1; j < q; j++) {
      const population = streamed[iq+j];
      newrho += population;
      const v = vectors[j];
      ${ Array.from({ length: d }, (_, i) => `u${i} += v[${i}]*population;`).join(' ') }
    }

    if (newrho <= 0) newrho = 0.01;

    rho[i] = newrho;
    ${ Array.from({ length: d }, (_, i) => `u${i} /= newrho;`).join(' ') }

    const tau = 3*viscosity + 0.5; // relaxation timescale
    const omega = 1 / tau;
    const invomega = 1 - omega;

    const usqr =  1 - 1.5 * (${ Array.from({ length: d }, (_, i) => `u${i}*u${i}`).join('+') });
    streamed[iq] = omega * weights[0] * newrho * usqr + invomega * streamed[iq];
    let dir;
    ${ expandVectorLoop(d, q, vectors, weights) }`;

  return new Function(
    'i', 'iq', 'viscosity', 'q', 'vectors', 'weights', 'rho', 'streamed', 'collided',
    body
  );
}

export default class LatticeBoltzmann {
  public rho: Float32Array;       // macroscopic density; cached for rendering
  private streamed: Float32Array; // microscopic densities along each lattice direction
  private collided: Float32Array;
  private size: number;           // total number of lattice points
  private dims: number[];
  private vectors: number[][];
  private weights: number[];
  private opp: number[];
  private d: number;
  private q: number;
  private _collide: Function;
  private _destination: Function;

  constructor(lattice: Lattice) {
    const dims = this.dims = lattice.dims;
    const d = this.d = dims.length;
    const weights = this.weights = lattice.weights;
    const q = this.q = weights.length;
    const size = this.size = dims.reduce((a, n) => a * n, 1);
    const vectors = this.vectors = lattice.vectors;
    this.streamed = createVectors(size, weights, q);
    this.collided = createVectors(size, weights, q);
    this.rho = new Float32Array(size);
    this.opp = lattice.opposites;
    this._collide = collide(d,q, vectors, weights);
    this._destination = destination(d, q, dims, size, vectors);
  }

  public step(viscosity: number, barriers: boolean[]) {
    const { size, _collide, _destination, vectors, weights, opp, q, rho, collided, streamed } = this;

    for (let i=0,iq=0; i<size; i++,iq+=q) {
      _collide(i, iq, viscosity, q, vectors, weights, rho, streamed, collided);
    }
    
    for (let i=0,iq=0; i<size; i++,iq+=q) {
      stream(i, iq, barriers, q, opp, _destination, collided, streamed);
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
    const usqr =  1 - 1.5 * dot(u, u, d);
    for (let j = 0; j < q; j++) {
      const dir = dot(vectors[j], u, d);
      streamed[iq+j] = weights[j] * rho * (usqr + 3 * dir + 4.5 * dir * dir);
    }
  }
}