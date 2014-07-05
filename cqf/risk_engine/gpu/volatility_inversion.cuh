enum kernel_t { serial, parallel, vega_parallel };

float call_price(const float s, const float r, const float v, const float t, const float k);
float call_vega(const float s, const float r, const float v, const float t, const float k);

void volatility_inversion(const float s, const float r, float *v, const float *t, const float *k, 
	const float *p, const float tol, const int iter, const int num, const kernel_t kernel);
