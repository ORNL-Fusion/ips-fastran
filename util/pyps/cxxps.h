/**
 * $Id: cxxps.h 855 2011-02-18 20:03:30Z pletzer $
 *
 */

#ifndef CXXPS_H
#define CXXPS_H

// local include 
#include "cpp/ccps.h"
#include "cpp/multiarray.h"

// std includes
#include <string>
#include <vector>
#include <map>

class PlasmaStateError {
 public:
  PlasmaStateError(const std::string& fN, int lN, int ier) :
    fileName(fN), lineNumber(lN), errorCode(ier) {};
  ~PlasmaStateError() {};
 private:
  std::string fileName;
  int lineNumber;
  int errorCode;
};

/**
 * A class that embodies all kinds of toroidal plasma data including 
 * geometry,  profiles, and hardware
 */
class PlasmaState {

public:

/**
 * Constructor
 */
    PlasmaState(const std::string& tag, const int debug);

/**
 * Destructor
 */
 virtual ~PlasmaState();

/** 
 * Return object handle 
 * @return handle
 */
 int* getHandle();

/**
 * Data setter
 * @param nm name of data
 * @param vl value(s) to set
 */
   void setData(const std::string& nm, const double vl);

/**
 * Data setter
 * @param nm name of data
 * @param vl value(s) to set
 */
   void setData(const std::string& nm, const int vl);

/**
 * Data setter
 * @param nm name of data
 * @param vl value(s) to set
 */
   void setData(const std::string& nm, const std::string& vl);

 /**
 * Data setter
 * @param nm name of data
 * @param vl value(s) to set
 */
   void setData(const std::string& nm, const std::vector<double>& vl);

/**
 * Data setter
 * @param nm name of data
 * @param vl value(s) to set
 */
   void setData(const std::string& nm, const std::vector<int>& vl);

/**
 * Data setter
 * @param nm name of data
 * @param vl value(s) to set
 */
    void setData(const std::string& nm, const std::vector<std::string>& vl);

/**
 * Data setter
 * @param nm name of data
 * @param vl value(s) to set
 */
    void setData(const std::string& nm, const MultiArray<double>& vl);

/**
 * Data setter
 * @param nm name of data
 * @param vl value(s) to set
 */
    void setData(const std::string& nm, const MultiArray<int>& vl);

/**
 * Data getter 
 * @param nm name of data
 * @return value(s)
 */
   template <typename T>
   T getData(const std::string& nm);
#ifdef SWIG
   %template(getDataInt)   getData<int>;
   %template(getDataDouble)   getData<double>;
   %template(getDataString)   getData<std::string>;
   %template(getDataDoubleVec)   getData< std::vector<double> >;
   %template(getDataStringVec)   getData< std::vector<std::string> >;
   %template(getDataIntVec)   getData< std::vector<int> >;
   %template(getDataMultiArrayDouble)   getData< MultiArray<double> >;
#endif

/** 
 * Get dimensions
 * @param nm name of data
 * @return dimensions
 */
   std::vector<size_t> getDims(const std::string& nm);

/**
 * Get underlying Fortran character size
 * @param nm name of data
 * @return size of CHARACTER declaration
 */
   size_t getCharSize(const std::string& nm);

/**
 * Get version Id
 * @return version_d
 */
   std::string getVersionId();

/**
 * Initial allocation of plasma state. Arrays will be filled with zeros.
 * Strings will be filled with blanks.
 */
   void alloc();

/** 
 * Set thermal species attributes
 * @param zAtom number of protons (atomic number)
 * @param zCharge 1 <= ionic charge <= zAtom
 * @param amu atomic mass unit
 */
   void setThermalSpecies(const int zAtom, const int zCharge, const int amu);

/** 
 * Set neutral beam species attributes
 * @param zAtom number of protons (atomic number)
 * @param zCharge 1 <= ionic charge <= zAtom
 * @param amu atomic mass unit
 */
   void setNeutralBeamSpecies(const int zAtom, const int zCharge, const int amu);

/** 
 * Set RF minority species attributes
 * @param zAtom number of protons (atomic number)
 * @param zCharge 1 <= ionic charge <= zAtom
 * @param amu atomic mass unit
 */
   void setRFMinoritySpecies(const int zAtom, const int zCharge, const int amu);

/** 
 * Set fusion species attributes
 * @param zAtom number of protons (atomic number)
 * @param zCharge 1 <= ionic charge <= zAtom
 * @param amu atomic mass unit
 */
   void setFusionSpecies(const int zAtom, const int zCharge, const int amu);

/** 
 * Set impurity species attributes
 * @param zAtom number of protons (atomic number)
 * @param zCharge 1 <= ionic charge <= zAtom
 * @param amu atomic mass unit
 */
   void setImpuritySpecies(const int zAtom, const int zCharge, const int amu);

/** 
 * Label, merge, and create neutrals
 */
   void finishSpecies();

/** 
 * Update MHD equilibrium stored in state object, from G-eqdsk file
 * 
 * @param g_filepath path to GEQDSK file
 * @param bdy_crat curvature ratio limit for boundary
 * @param kcur_option 0 or 1; 1=q(psi) computed from psirz, otherwise set from file
 * @param rho_curbrk rho value beyond which current density is interpolated to match total current
 */
   void updateEquilibrium(const std::string g_filepath, 
			  const double bdy_crat, 
			  const int kcur_option,
			  const double rho_curbrk);

/** 
 * Verify that profile of enclosed toroidal flux is consistent with rho coordinate
 *
 * @param reltol demanded relative accuracy (~0.02)
 * @param rhowk working grid (use internal grid if size is zero)
 * @param update_phit set to false if toroidal flux should not be updated
 * @param update_psi set to false if poloidal flux should not be updated
 * @param update_q set to false if safety factor should not be updated (update_psi and update_q cannot be both true)
 * @param check_phit set to true if toroidal flux should be checked
 * @param check_q set to true is safety factor should be checked
 * @return relerr_max max relative error in phi or q tests
 */
   double verifyMhdEq(const std::vector<double>& rhowk, const double reltol=0.02, 
	 const bool update_phit=false, 
         const bool update_psi=false, 
         const bool update_q=false, 
	 const bool check_phit=true, 
         const bool check_q=true);

/**
 * Derive flux surface averages and profiles from equilibrium geometry and fields
 * @param action specify action to be taken (eg "Everything")
 * @param rhowk working grid (use internal grid if size is zero)
 */
   void deriveMhdEq(const std::string& action, const std::vector<double>& rhowk);

/**
 * Read file containing machine description
 * @param filepath null terminated file path string (max 256 chars)
 * @param action erase plasma state if set to "NEW" or "INIT"
 * @param g_filepath null terminated file path string, to override box size and limiter contour (max 256 chars)
 * @param ierr error code (0=OK)
 */
   void readMachDescr(const std::string& filepath, const std::string& action,
		      const std::string& g_filepath="");

/**
 * Write commented machine description namelist from contents of state
 * @param filename null terminated file name (max 256 chars)
 */
   void writeMachDescr(const std::string& filename);

/**
 * Write geqdsk file based on contents of plasma state
 * @param filename null terminated file name (max 256 chars)
 */
   void writeGeqdsk(const std::string& filename);

/**
 * Update hash table; lock (or unlock) initialization sections of state
 * @param mdescr_flag if !=false => lock machine description section
 * @param sconfig_flag if !=false => lock shot description section
 * @param simInit_flag if !=false => lock simulation description section
 */
void updateHashCodes(const bool mdescr_flag, const bool sconfig_flag, 
           const bool simInit_flag);

/**
 * Read file containing shot configuration data
 * @param filepath file path string (max 256 chars)
 */
void readShotConfig(const std::string& filepath);

/**
 * Write file containing shot configuration data
 * @param filepath file path (max 256 chars)
 */
void writeShotConfig(const std::string& filepath);

/**
 * Write state to file
 * @param filepath  null terminated file path string (max 256 chars)
 */
 void store(const std::string& filepath);

/**
 * Summarize difference between two states, based on hash code comparison.
 * @param lbl1 label for 1st state
 * @param obj2 reference object
 * @param lbl2 label for 2nd state
 * @param sstr show all diffs for section, profile, dimension containg this string
 * @param icomp number of hash comparisons done (out)
 * @return idiff number of DIFFERENCES found
 */
 int hashDiff(const std::string& lbl1, 
	      PlasmaState& obj2, const std::string& lbl2,
	      const std::string& sstr, int* icomp);

/**
 * Read state object contents from file
 * @param filepath null terminated file path string (max 256 chars)
 * @param ierr error code (0=OK)
 */
 void getState(const std::string& filepath);

/** 
 * Interpolate 1d state element profile 
 * @param nm name of data
 * @param deriv 0 for value, 1 for d/dx, 2 for d^2/dx^2
 * @param xs target grid points
 * @return interpolated values
 */
 std::vector<double> interp1d(const std::string& nm, const int deriv, 
			      const std::vector<double>& xs);

/** 
 * Interpolate 1d state multi-element profile 
 * @param nm name of multi-data
 * @param deriv 0 for value, 1 for d/dx, 2 for d^2/dx^2
 * @param xs target grid points
 * @param index ID index in the multi-element profile
 * @return interpolated values
 */
 std::vector<double> interp1d(const std::string& nm, const int deriv, 
			      const std::vector<double>& xs, const size_t index);
/** 
 * Interpolate 2d state element profile f(x1, x2)
 * @param nm name of data
 * @param deriv1 0 for value, 1 for d/dx1, 2 for d^2/d[x1]^2
 * @param deriv2 0 for value, 1 for d/dx2, 2 for d^2/d[x2]^2
 * @param x1s target grid points
 * @param x2s target grid points (same size vector as above)
 * @return interpolated values
 */
 std::vector<double> interp2d(const std::string& nm, const int deriv1, 
			      const int deriv12,
			      const std::vector<double> x1s, 
			      const std::vector<double> x2s);

/** 
 * Interpolate 2d state multi-element profile f(x1, x2)
 * @param nm name of multi-data
 * @param deriv1 0 for value, 1 for d/dx1, 2 for d^2/d[x1]^2
 * @param deriv2 0 for value, 1 for d/dx2, 2 for d^2/d[x2]^2
 * @param x1s target grid points
 * @param x2s target grid points (same size vector as above)
 * @param index ID index in the multi-element profile
 * @return interpolated values
 */
 std::vector<double> interp2d(const std::string& nm, const int deriv1, 
			      const int deriv12,
			      const std::vector<double> x1s, 
			      const std::vector<double> x2s, 
			      const size_t index);

/**
 * Write file containing updates
 * @param filename null terminated file path string (max 256 chars)
 * @param hashflag != 0 to update I/O hash table
 */
 void writeUpdateFile(const std::string& filename, const int hashflag);

private:

/** opaque handle */
   int _h[12];

/** string -> string setter */
   std::map<std::string, 
     void (*)(int*, const char*, int*, size_t) > _setString;

/** string -> string getter */
   std::map<std::string, 
     void (*)(int*, char*, int*, size_t) >       _getString;

/** string -> string 1d setter */
   std::map<std::string, 
     void (*)(int*, const int*, const char*, int*, size_t) > _setString1D;

/** string -> string 1d getter */
   std::map<std::string, 
     void (*)(int*, const int*, char*, int*, size_t) >       _getString1D;

/** string -> int scalar setter */
   std::map<std::string, 
     void (*)(int*, const int*, int*) > _setInt;
/** string -> int scalar getter */
   std::map<std::string, 
     void (*)(int*, int*, int*) >       _getInt;

/** string -> int array 1d setter */
    std::map<std::string, 
      void (*)(int*, const int*, const int*, int*) > _setInt1D;

/** string -> int array 1d getter */
   std::map<std::string, 
     void (*)(int*, const int*, int*, int*) >       _getInt1D;

/** string -> int array 2d setter */
   std::map<std::string, 
     void (*)(int*, const int*, const int*, const int*, int*) > _setInt2D;
/** string -> int array 2d getter */
   std::map<std::string, 
     void (*)(int*, const int*, const int*, int*, int*) >       _getInt2D;

/** string -> int array 3d setter */
   std::map<std::string, 
     void (*)(int*, const int*, const int*, const int*, const int*, int*) > _setInt3D;
/** string -> int array 3d getter */
   std::map<std::string, 
     void (*)(int*, const int*, const int*, const int*, int*, int*) >       _getInt3D;

/** string -> double scalar setter */
   std::map<std::string, 
     void (*)(int*, const double*, int*) > _setDouble;

/** string -> double scalar getter */
   std::map<std::string, 
     void (*)(int*, double*, int*) >       _getDouble;

/** string -> double array 1d setter */
    std::map<std::string, 
      void (*)(int*, const int*, const double*, int*) > _setDouble1D;

/** string -> double array 1d getter */
   std::map<std::string, 
     void (*)(int*, const int*, double*, int*) >       _getDouble1D;

/** string -> double array 2d setter */
   std::map<std::string, 
     void (*)(int*, const int*, const int*, const double*, int*) > _setDouble2D;

/** string -> double array 2d getter */
   std::map<std::string, 
     void (*)(int*, const int*, const int*, double*, int*) >       _getDouble2D;

/** string -> double array 3d setter */
   std::map<std::string, 
     void (*)(int*, const int*, const int*, const int*, const double*, int*) > _setDouble3D;

/** string -> double array 3d getter */
   std::map<std::string, 
     void (*)(int*, const int*, const int*, const int*, double*, int*) >       _getDouble3D;

/** string -> double array 3d setter */
   std::map<std::string, 
     void (*)(int*, const int*, const int*, const int*, const int*, const double*, int*) > _setDouble4D;

/** string -> double array 3d getter */
   std::map<std::string, 
     void (*)(int*, const int*, const int*, const int*, const int*, double*, int*) >       _getDouble4D;

/** string -> array rank */
   std::map<std::string, 
     void (*)(int*, int*, int*) > _getNumDims;

/** string -> array dimensions */
   std::map<std::string, 
     void (*)(int*, const int*, int*, int*) > _getDims;

/** string -> get character size */
    std::map<std::string, void (*)(int*, int*, int*) > _getCharacterSize;

/** thermal species' atomic charges in MKS units */
    std::vector<double> _thAtomicCharges;
    
/** thermal species' ionic charges in MKS units */
    std::vector<double> _thIonicCharges;

/** thermal species' masses in kg */
    std::vector<double> _thMasses;

/** neutral beam species' atomic charges in MKS units */
    std::vector<double> _nbAtomicCharges;
    
/** neutral beam species' ionic charges in MKS units */
    std::vector<double> _nbIonicCharges;

/** neutral beam species' masses in kg */
    std::vector<double> _nbMasses;

/** RF minority species' atomic charges in MKS units */
    std::vector<double> _rfAtomicCharges;
    
/** RF minority species' ionic charges in MKS units */
    std::vector<double> _rfIonicCharges;

/** RF minority species' masses in kg */
    std::vector<double> _rfMasses;

/** fusion species' atomic charges in MKS units */
    std::vector<double> _fsAtomicCharges;
    
/** fusion species' ionic charges in MKS units */
    std::vector<double> _fsIonicCharges;

/** fusion species' masses in kg */
    std::vector<double> _fsMasses;

/** impurity species' atomic charges in MKS units */
    std::vector<double> _imAtomicCharges;
    
/** impurity species' ionic charges in MKS units */
    std::vector<double> _imIonicCharges;

/** impurity species' masses in kg */
    std::vector<double> _imMasses;
};

#endif // CXXPS_H
