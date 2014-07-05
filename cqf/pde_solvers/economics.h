#ifndef __ECONOMICS_H__
#define __ECONOMICS_H__

template<class T>
class equity_economics
{
    public :
        equity_economics(const T s, const T r, const T sigma)
            : _s(s), _r(r), _sigma(sigma) {  };
        
        T get_spot()            const { return _s;      }
        T get_interest_rate()   const { return _r;      }
        T get_volatility()      const { return _sigma;  }
        
    private :
        const T _s;
        const T _r;
        const T _sigma;
};

#endif
