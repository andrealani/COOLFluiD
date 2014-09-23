
#include "Environment/ObjectProvider.hh"

#include "StructMech/StructMech.hh"
#include "StructMech/MaterialLibKevlar.hh"

//////////////////////////////////////////////////////////////////////////////


namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MaterialLibKevlar,
                            MaterialPropertyLib,
                            StructMechModule,
                            1>
MaterialLibKevlarProvider("Kevlar");

//////////////////////////////////////////////////////////////////////////////

void MaterialLibKevlar::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal > ("YoungModulus", "Young Modulus");
  options.addConfigOption< CFreal > ("Poisson", "Poisson Coeficient");
  options.addConfigOption< CFreal > ("Density", "Density");
  options.addConfigOption< CFreal > ("ThermalExpansion", "Thermal Expansion Coef");

   options.addConfigOption< CFreal >("MaterialAngleX","Euler Angles in Deg between material coordinates and global coordinatesfor Anisotropic Materials");
   options.addConfigOption< CFreal >("MaterialAngleY","Euler Angles in Deg between material coordinates and global coordinatesfor Anisotropic Materials");
   options.addConfigOption< CFreal >("MaterialAngleZ","Euler Angles in Deg between material coordinates and global coordinatesfor Anisotropic Materials");

}

//////////////////////////////////////////////////////////////////////////////

MaterialLibKevlar::MaterialLibKevlar(const std::string& name) :
MaterialPropertyLib(name),
m_materialAngles(3),
m_T(6,6),
m_Tt(6,6)
{
  addConfigOptionsTo(this);

  // defaults are from
  m_young     = 1.;
  m_poisson   = 1.;
  m_density   = 1.;
  m_alpha     = 1.;
  m_isAnisotropic = true;

  ///@todo these should be vectors
  setParameter("YoungModulus", &m_young);
  setParameter("Poisson", &m_poisson);
  setParameter("Density", &m_density);
  setParameter("ThermalExpansion", &m_alpha);

  m_materialAngles = 0.;
  setParameter("MaterialAngleX",&m_materialAngles[0]);
  setParameter("MaterialAngleY",&m_materialAngles[1]);
  setParameter("MaterialAngleZ",&m_materialAngles[2]);

}

//////////////////////////////////////////////////////////////////////////////

MaterialLibKevlar::~MaterialLibKevlar()
{
}

//////////////////////////////////////////////////////////////////////////////

void MaterialLibKevlar::configure ( Config::ConfigArgs& args )
{
  MaterialPropertyLib::configure(args);

  m_C.resize(6,6);

  // Abaqus manual (taken from Kawabata and al. (1993))
  // E and G values are in N/m2
  CFreal E1 = 129600E6;
  CFreal E2 = 2490E6;
  CFreal E3 = 2490E6;

  CFreal nu12 = 0.31;
  CFreal nu23 = 0.0119;
  CFreal nu13 = 0.62;

  CFreal G12 = 2010E6;
  CFreal G23 = 924E6;
  CFreal G13 = 2010E6;

  // First set all values to zero
  m_C = 0;

  // Only the Compliance matrix is known
  // so invert it to get the stiffness matrix
  RealMatrix S(3,3);
  RealMatrix inverse(3,3);

  // Set non-zero values
  S(0,0)=1/E1;
  S(0,1)=-nu12/E1;
  S(0,2)=-nu13/E1;
  S(1,0)=S(0,1);
  S(1,1)=1/E2;
  S(1,2)=-nu23/E2;
  S(2,0)=S(0,2);
  S(2,1)=S(1,2);
  S(2,2)=1/E3;

  ///@todo change this, very bad (determinant is <<)
  m_inverter3.invert(S,inverse);
  for (CFuint i=0;i<3;++i){
    for (CFuint j=0;j<3;++j){
      m_C(i,j) = inverse(i,j);
    }
  }

  m_C(3,3)=G12;
  m_C(4,4)=G23;
  m_C(5,5)=G13;

  if(m_isAnisotropic){
    stiffnessMatrixTransform();
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibKevlar::computeYoungModulus()
{
  cf_assert(false);
  return m_young;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibKevlar::computePoissonCoef()
{
  cf_assert(false);
  return m_poisson;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibKevlar::computeThermalExpansionCoef()
{
  cf_assert(false);
  return m_alpha;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MaterialLibKevlar::computeDensity()
{
  return m_density;
}

//////////////////////////////////////////////////////////////////////////////

void MaterialLibKevlar::stiffnessMatrixTransform()
{

  /// Check the Inputed Euler Angles
  cf_assert(m_materialAngles.size() == 3);
  cf_assert(m_materialAngles[0] < 360.);
  cf_assert(m_materialAngles[1] < 360.);
  cf_assert(m_materialAngles[2] < 360.);

  /// Transform the angles into Radiants
  for (CFuint j=0; j < 3; ++j){
    m_materialAngles[j] *= MathTools::MathConsts::CFrealPi();
    m_materialAngles[j] /= 180.;
  }

  /// Transform the Anisotropic Stiffness Matrix into Global coordinates
  /// General Rotation Matrix given the Euler Angles
  // define temporary values
  CFreal c1 = cos(m_materialAngles[0]);
  CFreal c2 = cos(m_materialAngles[1]);
  CFreal c3 = cos(m_materialAngles[2]);
  CFreal s1 = sin(m_materialAngles[0]);
  CFreal s2 = sin(m_materialAngles[1]);
  CFreal s3 = sin(m_materialAngles[2]);

  // Define the vector rotation matrix components
  CFreal l1 = -s1*s3 + c1*c2*c3;
  CFreal l2 = c1*s3 + s1*c2*c3;
  CFreal l3 = -s2*c3;
  CFreal m1 = -s1*c3 - c1*c2*s3;
  CFreal m2 = c1*c3 - s1*c2*s3 ;
  CFreal m3 = s2*s3;
  CFreal n1 = c1*s2;
  CFreal n2 = s1*s2;
  CFreal n3 = c2;

  // Define the tensor rotation matrix
  ///@todo check order of rows and columns
  m_T(0,0) = l1*l1;
  m_T(0,1) = m1*m1;
  m_T(0,2) = n1*n1;
  m_T(0,3) = 2*l1*m1;
  m_T(0,4) = 2*m1*n1;
  m_T(0,5) = 2*n1*l1;

  m_T(1,0) = l2*l2;
  m_T(1,1) = m2*m2;
  m_T(1,2) = n2*n2;
  m_T(1,3) = 2*l2*m2;
  m_T(1,4) = 2*m2*n2;
  m_T(1,5) = 2*n2*l2;

  m_T(2,0) = l3*l3;
  m_T(2,1) = m3*m3;
  m_T(2,2) = n3*n3;
  m_T(2,3) = 2*l3*m3;
  m_T(2,4) = 2*m3*n3;
  m_T(2,5) = 2*n3*l3;

  m_T(3,0) = l1*l2;
  m_T(3,1) = m1*m2;
  m_T(3,2) = n1*n2;
  m_T(3,3) = l1*m2 + l2*m1;
  m_T(3,4) = m1*n2 + m2*n1;
  m_T(3,5) = n1*l2 + n2*l1;

  m_T(4,0) = l2*l3;
  m_T(4,1) = m2*m3;
  m_T(4,2) = n2*n3;
  m_T(4,3) = l2*m3 + l3*m2;
  m_T(4,4) = m2*n3 + m3*n2;
  m_T(4,5) = n2*l3 + n3*l2;

  m_T(5,0) = l3*l1;
  m_T(5,1) = m3*m1;
  m_T(5,2) = n3*n1;
  m_T(5,3) = l3*m1 + l1*m3;
  m_T(5,4) = m3*n1 + m1*n3;
  m_T(5,5) = n3*l1 + n1*l3;

  /// Compute the transposed matrix
  m_T.transpose(m_Tt);

  /// Modify the anisotropic stiffness matrix
  RealMatrix temp(6,6);
  temp = m_C * m_Tt;
  m_C = m_T * temp;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
