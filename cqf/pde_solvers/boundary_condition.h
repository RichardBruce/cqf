#ifndef __BOUNDARY_CONDITION_H__
#define __BOUNDARY_CONDITION_H__


template<class T>
class boundary_condition
{
    public :
        /* CTOR */
        boundary_condition(const T max_s, const T min_s)
         : _max_s(max_s), _min_s(min_s) {  };
         
        /* DTOR */
        virtual ~boundary_condition() {  };
         
        
        /* Access functions */
        T upper_boundary() const { return _max_s; }
        T lower_boundary() const { return _min_s; }
        
        
        /* Virtual apply */ 
        virtual void apply(T coeffs[3], const T s, const T r, const T sigma) const = 0;
        
    private :
        const T _max_s;
        const T _min_s;
};


template<class T>
class const_ds_boundary_condition : public boundary_condition<T>
{
    public :
        const_ds_boundary_condition(const T max_s, const T min_s)
         : boundary_condition<T>(max_s, min_s) {  };
         
        void apply(T coeffs[3], const T s, const T r, const T sigma) const
        {
            coeffs[0] = s * r;
            coeffs[1] = r - s * r;
            coeffs[2] = 0.0;
        }
};


template<class T>
class barrier_boundary_condition : public boundary_condition<T>
{
    public :
        barrier_boundary_condition(const T max_s, const T min_s)
         : boundary_condition<T>(max_s, min_s) {  };
         
        void apply(T coeffs[3], const T s, const T r, const T sigma) const
        {
            coeffs[0] = 0.0;
            coeffs[1] = 0.0;
            coeffs[2] = 0.0;
        }
};

#endif
