#pragma once


// -------------------------------------------------------------------------------------------- //
//models for calculation below-cloud and incloud scavenging coefficient, (1\s)
// -------------------------------------------------------------------------------------------- //


#include <math.h>
#include <stdio.h>
#include <string>
namespace nse
{
	template< typename T >
	class sc_coeff : public model_sc < T > {
	
	private:
	 
	 T ic_gr( const T pp_int, const T a_gr,const T b_gr) const; // Groell et al., 2014
	 T ic_ir( const T pp_int, const T a_ir,const T b_ir) const; // A. Querel, D. Quelo, Y. Roustan, A. Mathieu, 2021
	 T ic_na( const T pp_int, const T a_na,const T b_na) const; // Leadbetter et al., 2015
	 T ic_el( const T pp_int, const T a_el,const T b_el) const; // Querel et al, 2014
     T ic_en( const T pp_int, const T a_en,const T b_en) const; // Environ, 2018
	 T ic_sc( const T pp_int, const T a_sc,const T b_sc) const; // Scott, 1982
	 T ic_he(const T cl, const T fnuc ) const; // Hertel et al, 1995, Arnold et al, 2015, MRI, 2015, for 500m<Hi<2500m
	 T ic_ro(const T cl, const T deltaz, const T pp_int ) const; //	Roselle and Binkowski, 1999
	 T ic_sp(const T washout, const T t_cl) const; // B. Sportisse et al., 2006
	// T ic_ka(const T cl, const T fccn, const T fice ) const; // G. Katata, M. Chino
	 T fcl(const T pp_int) const;	
	 T ftwashout (const T LWC, const T deltaz, const T pp_int, const T ro_w) const;	
    // -------------------------------------------------------------------------------------------- //
     T bc_he( const T pp_int, const T a_he,const T b_he) const; // Arnold et al., 2015; Hertel et al., 1995
	 T bc_le( const T pp_int, const T a_le,const T b_le) const; // Leadbetter at al, 2015
	 T bc_mr( const T pp_int, const T a_mr,const T b_mr) const;// MRI, 2015
	 T bc_gr( const T pp_int, const T a_gl,const T b_gl) const;// Groell et all, 2014
	
	public:
		// data (public):
		T pp_int; // precipitation intensity (mm/hr)
        T fnuc; // the fraction of active aerosols, may be given by the constant 0.9
	    T fccn; //
		T fice; //
		T deltaz; //cloud thickness (m)
	    T LWC; //cloud liqued water
		T ro_w; // water density
	    T cl;
	    T twashout;
	private:

		// fixed coefficient aI^b
		static const T a_gr;
		static const T a_ir;
		static const T a_na;
		static const T a_el;
		static const T a_en;
		static const T a_cs;
		static const T b_gr;
		static const T b_ir;
		static const T b_na;
		static const T b_el;
		static const T b_en;
		static const T b_gr;
	    static const T t_cl;
	    static const T a_gl;
		static const T b_gl;
	    static const T a_mr;
		static const T b_mr;
		static const T a_he;
		static const T b_he;
	    static const T a_le;
		static const T b_le;
	
	
	};
}
// -------------------------------------------------------------------------------------------- //

//fixed coefficient aI^b

template< typename T >
const  T nse::sc_coeff< T >::a_gr = (T)0.00005;
template< typename T >
const  T nse::sc_coeff< T >::a_ir = (T)0.0005;
template< typename T >
const  T nse::sc_coeff< T >::a_na = (T)0.000336;
template< typename T >
const  T nse::sc_coeff< T >::a_el = (T)0.000397;
template< typename T >
const  T nse::sc_coeff< T >::a_en = (T)0.00042;
template< typename T >
const  T nse::sc_coeff< T >::a_sc = (T)0.00035;
template< typename T >
const  T nse::sc_coeff< T >::b_gr = (T)1;
template< typename T >
const  T nse::sc_coeff< T >::b_ir = (T)0.64;
template< typename T >
const  T nse::sc_coeff< T >::b_na = (T)0.79;
template< typename T >
const  T nse::sc_coeff< T >::b_el = (T)0.31;
template< typename T >
const  T nse::sc_coeff< T >::b_en = (T)0.79;
template< typename T >
const  T nse::sc_coeff< T >::b_sc = (T)0.78;
template< typename T >
const  T nse::sc_coeff< T >::t_cl = (T)3600;
template< typename T >
const  T nse::sc_coeff< T >::b_he = (T)0.8;
template< typename T >
const  T nse::sc_coeff< T >::b_le = (T)0.79;
template< typename T >
const  T nse::sc_coeff< T >::b_mr = (T)0.75;
template< typename T >
const  T nse::sc_coeff< T >::b_gl = (T)1;
template< typename T >
const  T nse::sc_coeff< T >::a_he = (T)0.0004;
template< typename T >
const  T nse::sc_coeff< T >::a_le = (T)0.000084;
template< typename T >
const  T nse::sc_coeff< T >::a_mr = (T)0.000029;
template< typename T >
const  T nse::sc_coeff< T >::a_gl = (T)0.00005;
// -------------------------------------------------------------------------------------------- //
//  functions for models
// -------------------------------------------------------------------------------------------- //
	
template< typename T >
  T  nse::surface_snow< T >::ic_gr(  const T pp_int, const T a_gr,const T b_gr ) const 
	{
	  return a_gr*pow(pp_int,b_gr);
	}

template< typename T >
  T  nse::surface_snow< T >::ic_ir(  const T pp_int, const T a_ir,const T b_ir ) const 
	{
	  return a_ir*pow(pp_int,b_ir);
	}
	
	template< typename T >
  T  nse::surface_snow< T >::ic_na(  const T pp_int, const T a_na,const T b_na ) const 
	{
	  return a_na*pow(pp_int,b_na);
	}
	
	template< typename T >
  T  nse::surface_snow< T >::ic_el(  const T pp_int, const T a_el,const T b_el ) const 
	{
	  return a_el*pow(pp_int,b_el);
	}
	
	template< typename T >
  T  nse::surface_snow< T >::ic_en(  const T pp_int, const T a_en,const T b_en ) const 
	{
	  return a_en*pow(pp_int,b_en);
	}
	
	template< typename T >
  T  nse::surface_snow< T >::ic_sc(  const T pp_int, const T a_sc,const T b_sc ) const 
	{
	  return a_sc*pow(pp_int,b_sc);
	}



template< typename T >
  T  nse::surface_snow< T >::fcl(const T pp_int ) const 
	{
	  return (T)0.0000002*pow(pp_int,(T)0.36);
	}
		
template< typename T >
  T  nse::surface_snow< T >::ic_he(  const T cl, const T fnuc ) const 
	{
	  return fnuc/((T)3600000*cl);
	}
	
	template< typename T >
  T  nse::surface_snow< T >::ic_ro(  const T cl, const T deltaz, const T pp_int) const 
	{
	  return (T)1/(T)3600)*(1-pow(exp, ((T)-1000*deltaz*(cl/pp_int))));
	}
	
	template< typename T >
T  nse::surface_snow< T >:: ftwashout (const T LWC, const T deltaz, const T pp_int, const T ro_w)
		{
	  return (LWC*deltaz)/(ro_w*pp_int);
	}	
			
		template< typename T >
  T  nse::surface_snow< T >::ic_sp(const T twashout, const T t_cl) const 
	{
	  return ((T)1-pow(exp,(-t_cl/twashout))/t_cl;
	}
	
	
	template< typename T >
  T  nse::surface_snow< T >::bc_he(  const T pp_int, const T a_he,const T b_he ) const 
	{
	  return a_he*pow(pp_int,b_he);
	}
	
		template< typename T >
  T  nse::surface_snow< T >::bc_le(  const T pp_int, const T a_le,const T b_le ) const 
	{
	  return a_le*pow(pp_int,b_le);
	}
	
		template< typename T >
  T  nse::surface_snow< T >::bc_mr(  const T pp_int, const T a_mr,const T b_mr ) const 
	{
	  return a_mr*pow(pp_int,b_mr);
	}
	
		template< typename T >
  T  nse::surface_snow< T >::bc_gl(  const T pp_int, const T a_gl,const T b_gl ) const 
	{
	  return a_gl*pow(pp_int,b_gl);
	}
	
	
	
// -------------------------------------------------------------------------------------------- //
