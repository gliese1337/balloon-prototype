uniform sampler2D streamed[19]; // single channel
uniform vec3 c[19];

void main () {
  float rho = texture2D(streamed[0], gl_FragCoord).r;
  vec3 u = vec3(0.0); // macroscopic velocity
  for (int j = 1; j < 19; j++) {
    float v = texture2D(streamed[j], gl_FragCoord);
    rho += v;
    u += v*c[j];
  }

  // hack for stability
  rho = max(rho, 0.01);

  gl_FragColor = vec4(u / rho, rho);
}