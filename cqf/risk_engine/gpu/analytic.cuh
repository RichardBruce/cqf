__device__ float call_price(const float s, const float r, const float v, const float t, const float k);
__device__ float call_vega(const float s, const float r, const float v, const float t, const float k);

__device__ float put_price(const float s, const float r, const float v, const float t, const float k);
__device__ float put_vega(const float s, const float r, const float v, const float t, const float k);
