#include "XYEnsemble.h"
#include "XYSIS.h"
#include <vector>
//#include "debug_new.h"

using namespace spatialaggregate;
using std::vector;

//----------------------------------------------------------------------------
CXYEnsemble::CXYEnsemble()
{
}

CXYEnsemble::CXYEnsemble(
	char* cOutPath,
	float fPersistenLength,
	float fCollisionLength,
	float fPackingDensity,
	float fNucleusSphereDiameter,
	int iNumNodes,
	int iNumSamplePoints,
	int NumberofChains,
	vector <int> ChainLengths,
	char* cStartEndFile,
	char* cContIndFile)
{
	m_cOutPath = new char[1024];
	SetOutPath(cOutPath);
	SetPersistenceLength( fPersistenLength );
	SetPackingDensity(fPackingDensity);
	SetNucleusSphereDiameter(fNucleusSphereDiameter);
	SetNumChains(NumberofChains);
	SetChainLengths(ChainLengths);

	m_iNumSamplePoints = iNumSamplePoints;
	if (strcmp(cStartEndFile, "") == 0){
		// coarse version
		// given number of nodes
		SetNumNodes(iNumNodes);
		// set each segment length equal to persistence length
		SetSegLengths();
    SetContIndex();
	} else {
		// fine version
		// read start end position file
		// dynamic setting each segment length determined by mass density
		SetSegLengths(cStartEndFile, "AVG");
		SetNumNodes(m_vfSegLength.size()); // set node numbers
    SetContIndex(cContIndFile);
	}
	SetCollisionLength( fCollisionLength );

    mNum =10;   //number of the chains

	//---------------------------------------- GAMZE 10/2 -------------------------------------




	for (int i = 0; i < mNum; i++)
	{

	     vector <tree <CXYVector <float> >* > v;
         for (int j=0; j<m_fNumofChains; j++)
         {

             tree <CXYVector<float> >* TreeTemp= new tree< CXYVector <float> >  () ;
             v.push_back(TreeTemp);
            // delete TreeTemp;
         }
         m_pTChain.push_back(v);
        // v.clear();

	}







	//--------------------------------------------------------------------------------------------------

	m_pMSamplesOrg = new CXYMatrix<float>;
	SetSamplesOrg();

  m_iNumMiddleEndPoints = max(int(m_fPersistenceLength / m_fCollisionLength)-1 ,1);

  cout << "persistence length = " << m_fPersistenceLength << endl;
  cout << "collision length = " << m_fCollisionLength << endl;
  cout << "number nodes = " << m_iNumNodes << endl;
  cout << "m_iNumMiddleEndPoints = " << m_iNumMiddleEndPoints <<endl;


  // Construct octree
  m_center = Eigen::Matrix< float, 4, 1 > (0.0f,0.0f,0.0f,1);
  m_minimumVolumeSize = m_fCollisionLength/2;
  m_dr = m_fCollisionLength*2;
//  m_dr = m_fNucleusSphereDiameter/2;

//----------------------------------------------GAMZE 10/2--------------------------------------------------------------------------------------------------------------

    for (int i = 0; i < mNum; i++)
	{
        boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  =  boost::make_shared< OcTreeNodeAllocator< float , int > >();
        m_pOctree.push_back(new OcTree<float,int>(m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator) );

	}

	//logweight per sample

	  for (int i=0 ; i<mNum; i++)
        {
            //cout<<"index is "<<i<<endl;
            ms_logweight.push_back(FloatIntPair(0,i));

            //cout << "logweight is"<<m_LogWeight[i].first <<endl;
        }

        //logweight per chain

//  for (int i=0 ; i<m_fNumofChains; i++)
//        {
//            //cout<<"index is "<<i<<endl;
//            m_LogWeight.push_back(FloatIntPair(0,i));
//
//            //cout << "logweight is"<<m_LogWeight[i].first <<endl;
//        }

m_maxDepth = ceil(m_pOctree[0]->depthForVolumeSize(m_minimumVolumeSize));
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  cout << "minimumVolumeSize =" << m_minimumVolumeSize << endl;
  cout << "maxDistance =" << m_fNucleusSphereDiameter << endl;
  cout << "maxDepth = " <<  m_maxDepth << endl;


}

//----------------------------------------------------------------------------
CXYEnsemble::~CXYEnsemble(void)
{


    for (int i=0 ; i<mNum; i++)
    {
        delete  m_pOctree[i];
   //     delete  m_pTChain1[i];

    }


//
 for (int i=0 ; i<mNum; i++)
    {
        for (int j=0; j<m_fNumofChains; j++)
        {
            delete m_pTChain[i][j];
        }
            //delete  v[i];

    }

	delete m_cOutPath;
	delete m_pMSamplesOrg;
	//m_pOctree.empty();
    //delete  m_fChainLengths;//.empty();
//    m_pTChain.empty();
  //  m_LogWeight.empty();
  //  ms_logweight.empty();
}

//----------------------------------------------------------------------------
void CXYEnsemble::SetPersistenceLength(float fPL)
{
	m_fPersistenceLength = fPL;
}
//----------------------------------------------------------------------------
float CXYEnsemble::GetPersistenceLength(void)
{
	return m_fPersistenceLength;
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetNumChains(int NofC)
{
	m_fNumofChains = NofC;
}
//----------------------------------------------------------------------------
int CXYEnsemble::GetNumChains(void)
{
	return m_fNumofChains;
}
//------------------------------------------------------------------------------

void CXYEnsemble::SetChainLengths(vector <int> fCLs)
{


    int NofC = GetNumChains();

    for (int i=0; i<NofC; i++)
    {
        m_fChainLengths.push_back(fCLs[i]);
    }

}
//----------------------------------------------------------------------------
vector <int> CXYEnsemble::GetChainLengths(void)
{
	return m_fChainLengths;
}
//------------------------------------------------------------------------------
void CXYEnsemble::SetCollisionLength(float fCollision)
{
	vector<float>& vfSegment = GetSegLengths();
	float min_diameter = *(min_element(vfSegment.begin(), vfSegment.end()));
	m_fCollisionLength = min(fCollision,min_diameter);
}
//----------------------------------------------------------------------------
float CXYEnsemble::GetCollisionLength(void)
{
	return m_fCollisionLength;
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetPackingDensity(float fPackingDensity)
{
	m_fPackingDensity = fPackingDensity;
}
//----------------------------------------------------------------------------
float CXYEnsemble::GetPackingDensity(void)
{
	return m_fPackingDensity;
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetNucleusSphereDiameter(float fNucleusSphereDiameter)
{
	m_fNucleusSphereDiameter = fNucleusSphereDiameter;
}
//----------------------------------------------------------------------------
float CXYEnsemble::GetNucleusSphereDiameter(void)
{
	return m_fNucleusSphereDiameter;
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetNumNodes(int iNumNodes){
	m_iNumNodes = iNumNodes;
}
//----------------------------------------------------------------------------
int CXYEnsemble::GetNumNodes(void){
	return m_iNumNodes;
}

//----------------------------------------------------------------------------
// Generate sphere sample points with radius persistence length.
void CXYEnsemble::SetSamplesOrg(void)
{
	CXYSO3Sequence sO3sequence(m_iNumSamplePoints);
	sO3sequence.SetSO3Sequence();
	(*m_pMSamplesOrg) = sO3sequence.GetSO3Sequence();
}
//----------------------------------------------------------------------------
CXYMatrix<float> CXYEnsemble::GetSamplesOrg()
{
	return (*m_pMSamplesOrg);
}

//----------------------------------------------------------------------------
//CXYMatrix<float> CXYEnsemble::NewXYZ(CXYVector<float> kV_0, CXYVector<float> kV_1, CXYVector<float> kV_2)
//{
//	CXYMatrix<float> kMXYZ(3,3);
//	CXYVector<float> kV_Z = kV_2 - kV_1;
//	CXYVector<float> kV_01 = kV_1 - kV_0;
//	kV_01.Normalize();
//	if (fabs(kV_Z.Dot(kV_01)) < CXYMath<float>::ZERO_TOLERANCE)
//	{ // two line parallel
//		kMXYZ = NewXYZ(kV_1,kV_2);
//	} else {
//		float a = kV_Z[0],  b = kV_Z[1],  c = kV_Z[2];
//		float x = kV_01[0], y = kV_01[1], z = kV_01[2];
//
//		float fX[3] = {	y*c - z*b, z*a - x*c, x*b - y*a };
//		CXYVector<float> kV_X(3,fX);
//		kV_X.Normalize();
//
//		x = kV_Z[0], y = kV_Z[1], z = kV_Z[2];
//		a = kV_X[0], b = kV_X[1], c = kV_X[2];
//		float fY[3] = {	y*c - z*b, z*a - x*c, x*b - y*a };
//		CXYVector<float> kV_Y(3,fY);
//		kV_Y.Normalize();
//		kMXYZ.SetRow(0, kV_X);
//		kMXYZ.SetRow(1, kV_Y);
//		kMXYZ.SetRow(2, kV_Z);
//	}
//
//	return kMXYZ;
//}
////----------------------------------------------------------------------------
//CXYMatrix<float> CXYEnsemble::NewXYZ(CXYVector<float> kV_1, CXYVector<float> kV_2)
//{
//	CXYMatrix<float> kMXYZ(3,3);
//	CXYVector<float> kV_Z = kV_2 - kV_1;
//	kV_Z.Normalize();
//	float a = kV_Z[0], b= kV_Z[1], c = kV_Z[2];
//
//	float x, y, z;
//	if (c > CXYMath<float>::ZERO_TOLERANCE )
//	{
//		x = 1; y = 0; z = -(a*x+b*y)/c;
//	} else if (b > CXYMath<float>::ZERO_TOLERANCE )
//	{
//		x = 0; z = 1; y = -(a*x+c*z)/b;
//	} else {
//		y = 1; z = 0; x = -(b*y+c*z)/a;
//	}
//	float fX[3] = {x,y,z};
//	CXYVector<float> kV_X(3,fX);
//	float fY[3] = {b*z-c*y, c*x-a*z, a*y-b*x};
//	CXYVector<float> kV_Y(3,fY);
//
//	kMXYZ.SetRow(0, kV_X);
//	kMXYZ.SetRow(1, kV_Y);
//	kMXYZ.SetRow(2, kV_Z);
//
//
//	return kMXYZ;
//}

//----------------------------------------------------------------------------
//CXYMatrix<float> CXYEnsemble::GetRotMatrix(CXYMatrix<float> &rkMXYZ)
//{
//	CXYVector<float> kV_X = rkMXYZ.GetRow(0);
//	CXYVector<float> kV_Y = rkMXYZ.GetRow(1);
//	CXYVector<float> kV_Z = rkMXYZ.GetRow(2);
//
//	CXYMatrix<float> kM_A(3,3);
//	float alpha, beta, gamma;
//	float Z3diff = fabs(1-kV_Z[2]*kV_Z[2]);
//	if ( Z3diff < CXYMath<float>::ZERO_TOLERANCE)
//	{
//		alpha = acos( min(float(1.0), max(float(-1.0), kV_X[0])));
//		beta  = 0;
//		gamma = 0;
//	} else {		// http://www.macosxguru.net/article.php?story=20040210124637626
//		float Denorm = sqrt(Z3diff);
//		alpha = acos(min(float(1.0), max(float(-1.0), -kV_Z[1]/Denorm)));
//		beta  = acos(min(float(1.0), max(float(-1.0), -kV_Z[2])));
//		gamma = acos(min(float(1.0), max(float(-1.0), -kV_Y[2]/Denorm)));
//	}
//
//	kM_A[0][0] = cos(gamma) * cos(alpha) - cos(beta) * sin(alpha) * sin(gamma);
//	kM_A[0][1] = cos(gamma) * sin(alpha) + cos(beta) * cos(alpha) * sin(gamma);
//	kM_A[0][2] = sin(gamma) * sin(beta);
//	kM_A[1][0] = -sin(gamma) * cos(alpha) - cos(beta) * sin(alpha) * cos(gamma);
//	kM_A[1][1] = -sin(gamma) * sin(alpha) + cos(beta) * cos(alpha) * cos(gamma);
//	kM_A[1][2] = cos(gamma) * sin(beta);
//	kM_A[2][0] = sin(beta) * sin(alpha);
//	kM_A[2][1] = -sin(beta) * cos(alpha);
//	kM_A[2][2] = cos(beta);
//
//	kM_A.GetInverse(kM_A);
//
//
//	return kM_A;
//}

//----------------------------------------------------------------------------
void CXYEnsemble::InitializeChain(int m_samples)
{
  // renew octree


  //---------------------------------------------------GAMZE 10/2-----------------------------------------------

//        for ( int i = 0; i < m; i++)
//        {
//           // cout<<"index is "<<i<<endl;
//            m_LogWeight[i].first=0;
//
//           // cout << "logweight is"<<m_LogWeight[i].first <<endl;
//        }

        for ( int i = 0; i < m_samples; i++)
        {
          //  cout<<"index is "<<i<<endl;
            ms_logweight[i].first=0;

          //  cout << "logweight is"<<m_LogWeight[i].first <<endl;
        }


//if octree is not null, then delete it ------------------
  for ( int i = 0; i < m_pOctree.size(); i++ )
  {

       if (m_pOctree[i] -> root_) {
            delete m_pOctree[i];

       }
  }
m_pOctree.clear();
  //---------------------------------------------------------------------------------------------------------------------


  OcTree<float,int> * p_octree;



//----------------------here I am-------------------------


  vector < vector < tree < CXYVector <float> >*> > ptr = GetTree();



	// erase all nodes if not empty



for (int j=0; j < m_samples; j++)
{
  for (int i = 0; i < m_fNumofChains; i++)
  {
      //cout<<"initialize chain"<<endl;
      if ( !ptr[j][i]-> empty() )
      {
          ptr[j][i]-> clear();
      }
  }
}




vector < vector  <  tree < CXYVector <float> >::iterator> >  top, root;



for (int j=0; j<m_samples; j++)
{
    vector  <  tree < CXYVector <float> >::iterator> t1,r1;

    for (int i = 0; i < m_fNumofChains; i++)
    {
        tree < CXYVector <float> >::iterator TempI;
        t1.push_back(TempI);
        r1.push_back(TempI);
    }
    top.push_back(t1);
    root.push_back(r1);
}

for (int i=0; i<m_samples; i++)
{
    for (int j=0; j<m_fNumofChains; j++)
    {
        top[i][j]= ptr[i][j]-> begin();
    }
}




CXYVector<float>  kV_startpoint;
//cout << "Size = " << m_pOctree.size() << endl;
for (int j=0; j < m_samples; j++)
{
    boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  =  boost::make_shared< OcTreeNodeAllocator< float , int > >();
    m_pOctree.push_back(new OcTree<float,int>(m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator) );

    for ( int i = 0; i < m_fNumofChains; i++ )
    {
        kV_startpoint = RndSetStartPoint();
        root[j][i] = ptr[j][i]-> insert (top[j][i], kV_startpoint);

    // boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  =  boost::make_shared< OcTreeNodeAllocator< float , int > >();
    // p_octree = new OcTree<float,int> (m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator);
    // m_pOctree.push_back(p_octree);
        m_pOctree[j] ->addPoint(kV_startpoint[0],kV_startpoint[1],kV_startpoint[2],1,m_maxDepth);

     //cout<<"adress is "<<m_pOctree[i]<<endl;



   // m_pOctree.push_back(p_octree);
 //cout<<"address"<<*(m_pOctree[i])<<endl;
    }
}


cout<<"initialize chain has ended"<<endl;

}
//----------------------------------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------GAMZE 10/2-------------------------------------------
vector < vector < tree< CXYVector<float> >* >  > CXYEnsemble::GetTree()
{
    //cout<<"get tree"<<endl;
	return m_pTChain;
}
//------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------
CXYVector<float> CXYEnsemble::RndSetStartPoint(void)
{
	float fdiameter = GetNucleusSphereDiameter();
	float fradius = fdiameter/2;
	MTRand mtrand;

	// Check if the random point located in the neuclus sphere.
	// When one point satisfy the condition, we accept it.
	float fCoord[3];
	while (1) {
		fCoord[0] = mtrand.randExc(fradius);
		fCoord[1] = mtrand.randExc(fradius);
		fCoord[2] = mtrand.randExc(fradius);
		if (IsInsideSphere(fCoord)){
			break;
		}
	}
//	cout << x << ";" << y << ";" << z<<";" <<endl;
//	cout << x*x+y*y+z*z << endl;
//	cout << fradius2 << endl;
	CXYVector<float> kV (3,fCoord);
	return kV;
}


//----------------------------------------------------------------------------
bool CXYEnsemble::IsInsideSphere(CXYVector<float> kV_point)
{

    //cout<< "girdik mi ki buraya"<<endl;
	float fradius = GetNucleusSphereDiameter()/2;
	float fradius2 = fradius*fradius;
	bool flag;

	if (kV_point.SquaredLength() < fradius2) {
	  //   cout<<"girdik mi buraya"<<endl;

		flag = true;

	}else {
		flag = false;
	}

    //cout<<"ciktik burdan"<<endl;
   // cout<<"flag is "<<flag<<endl;
	return flag;
	//cout<<"ciktik burdan"<<endl;


}
//----------------------------------------------------------------------------
bool CXYEnsemble::IsInsideSphere(float* fCoord)
{
//    cout<<"girdik mi buraya"<<endl;
	float fradius = GetNucleusSphereDiameter()/2;
	float fradius2 = fradius*fradius;
	bool flag;
	if (fCoord[0]*fCoord[0] + fCoord[1]*fCoord[1] + fCoord[2]*fCoord[2]
			< fradius2) {
		flag = true;
	}else {
		flag = false;
	}
//	cout<<"ciktik mi burdan"<<endl;
	return flag;
}

//---------------------------------------GAMZE 10/11-------------------------------------
bool CXYEnsemble::IsCollision(CXYVector<float> kV_point, int j)
{

 // cout<<"WTF 3"<<endl;
  OcTree<float,int> * p_octree;

  float x_ = kV_point[0];
  float y_ = kV_point[1];
  float z_ = kV_point[2];
  //cout<<x_<<" "<<y_<<" "<<z_<<endl;

  float min_x_ = max(x_ - m_dr, -m_fNucleusSphereDiameter);
  float min_y_ = max(y_ - m_dr, -m_fNucleusSphereDiameter);
  float min_z_ = max(z_ - m_dr, -m_fNucleusSphereDiameter);
  Eigen::Matrix<float,4,1> minPosition_ (min_x_,min_y_,min_z_,1);
  Eigen::Matrix<float,4,1> maxPosition_ (x_+m_dr,y_+m_dr,z_+m_dr,1);
   vector< OcTreeNode< float, int >* > nodes_;


//cout<<"m_pOctree is"<<m_pOctree[j]<<endl;
//boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  =  boost::make_shared< OcTreeNodeAllocator< float , int > >();
 //    p_octree = new OcTree<float,int> (m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator);
 //    m_pOctree.push_back(p_octree);

   m_pOctree[j]->getAllNodesInVolumeOnDepth(nodes_,minPosition_,maxPosition_,m_maxDepth,true);

  // delete p_octree;

 //   cout<<"this is working"<<endl;


  if (nodes_.size() > 0) {
    bool bFlag = false;
    for (int i = 0; i< (int)nodes_.size(); i++) {
      Eigen::Matrix<float,4,1> point_ = (nodes_[i])->getPosition();
      float point_x_ = point_[0];
      float point_y_ = point_[1];
      float point_z_ = point_[2];
      float dist_ = sqrt((point_x_ - x_)*(point_x_ - x_) + (point_y_ - y_ )*(point_y_ - y_) + (point_z_ - z_)*(point_z_ - z_));
      if (dist_ < m_fCollisionLength) {
        bFlag = true;
        break;
//        cout << x_ << " " << y_ << " " << z_ <<endl;
//        cout << point_x_ << " " << point_y_ << " " << point_z_ <<endl;
//        cout << dist_ << endl;
//        cout << "collision happend";
      }
    }
    return bFlag;
	//return false;
  }else {
    return false;
  }
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------
//--------------------------------GAMZE 10/3------------------------------------------
bool CXYEnsemble::GrowmChain(int k, int f, int m_samples)
{

    //int f=0;

     bool Flag=true;
//-----------initialize m chain if it's already not initialized---------------------------------------------------------

    cout<<"grow m chain is starting"<<endl;
    if ( f == 0 )
    {


        cout<<"initialize chain is started"<<endl;
        InitializeChain(m_samples);
        cout<<"initialize chain is done"<<endl;

    }

//--------------------------------------------------------------------------------------------



	MTRand mtrand;

	int iNumNodes = GetNumNodes();
	vector <int> ChainLenghts = GetChainLengths();


//----------------------------vector of matrices---------------------------------------------
	vector <vector<tree< CXYVector<float> >*> > ptr = GetTree();


//---------------------------------------------------------------------------------------------





//-----------------------------initialize the weights----------------------------------



//---------------------------------------------------------------------------------------------
//  iNumNodes = 3;

//-------------------grow m chain-------------------------------------------------------
//-------------------each chain has k number of nodes------------------------
//two scenario: 1) if chain is already grown means f>0----------------------------------
//--------------------2) we just started to grow means f=0 ------------------------------------



    int a=f;
    int numnodes=(f*k) +2*k;
    cout<<"numnodes is "<<numnodes<<endl;


    int j=0;
    vector <int> pro_samp1;
    vector <int> pro_samp;

   for (int i=0; i<m_samples; i++)
   {
	pro_samp1.push_back(-1);
	pro_samp.push_back(-1);
   }

    int ind=0;
    int ind1=0;

    for (int i=0; i<m_samples; i++)
    {
     //   cout<<"this is happenning"<<endl;
        while(j<m_fNumofChains)
        {
            if (numnodes<ChainLenghts[j])
            {

                if (GrowOneChain(i,j,numnodes,a))
                {
                    j=j+1;
                    counter=0;
                }
                else
                {
                    j=j;
                    counter=counter+1;
                    if (counter >=500)
                    {
                        indicator=1;
                       // Flag = true;
			//ms_logweight[i].first=0;
                	    cout<<"counter is "<<counter<<endl;
    	                pro_samp1[ind1]=i;
        	            counter=0;
                        j=j+1;
                        ind1 = ind1+1;

                    }
                   // cout<<"counter is "<<counter<<endl;
                }
               // counter = 0;

            }
            else if (numnodes>=ChainLenghts[j] && (numnodes-k)<ChainLenghts[j])
            {
                if (GrowOneChain(i,j,ChainLenghts[j],a))
                {
                    j=j+1;
                    counter=0;
                }
                 else
                {
                    j=j;
                     counter=counter+1;
                    //cout<<"counter is "<<counter<<endl;
                    if (counter >=500)
                    {
		
                        cout<<"counter is "<<counter<<endl;
                        indicator = 1;
                       // Flag = true;
                        pro_samp[ind]=i;
		//	ms_logweight[i].first=0;
                        counter=0;
                        j=j+1;
                        ind=ind+1;
			exit(0);
                    }
                }

            }


            else
            {
                j=j+1;
            }
        }
        j=0;
    }

//cout<<"is it out of loop"<<endl;

  for (int i=0; i<ind; i++)
  {
      ms_logweight[pro_samp[i]].first = -5000;
  }

   for (int i=0; i<ind1; i++)
  {
      ms_logweight[pro_samp1[i]].first = -5000;
  }


 //cout<<"here I am"<<endl;
  return Flag;

}

bool CXYEnsemble::GrowOneChain(int m_samples, int j, int numnodes, int a)

{

        bool Flag=true;

        float sum = 0;

       vector <vector< tree< CXYVector <float> >* > > ptr = GetTree();

       MTRand mtrand;

//	int iNumNodes = GetNumNodes();

     //   counter = 0;
        for ( int i = 1+ a; i < numnodes; i++ )
        {
            //cout<<"is this happening"<<endl;
            tree< CXYVector<float> >::pre_order_iterator node;
            node = ptr[m_samples][j] -> end();
            node --;




            CXYMatrix < float > kMSamplesPoints = GetSamplesOrg();
            kMSamplesPoints *= GetSegLength(i-1);



            vector<int> GoodPointInd;
            vector<int> NoCollisionPointInd;

       //      cout<<"node is "<<(*node)<<endl;

            GetGoodPoints( (*node), kMSamplesPoints, GoodPointInd, NoCollisionPointInd,m_samples);

         //   cout<<" get good points HAPPENED"<<endl;
            int GoodPointSize = GoodPointInd.size();
            int NoCollisionPointSize = NoCollisionPointInd.size();
         //   cout<<"good point size is "<<GoodPointSize<<endl;
         //   cout<<"no collision point size is "<<NoCollisionPointSize<<endl;
            if (GoodPointSize == 0 || NoCollisionPointSize == 0) {
                Flag = false;
             //   counter = counter + 1;
                break;
            }




          // cout<<"logweight is"<<m_LogWeight[j].first<<endl;
     //       m_LogWeight[j].first += log(NoCollisionPointSize);
      //      m_vLogWeight.push_back(m_LogWeight[j].first);
           // sum = sum + m_LogWeight[j].first;
            ms_logweight[m_samples].first += log(NoCollisionPointSize);
         //   cout<<"logweight is "<<ms_logweight[m_samples].first<<" index is "<<ms_logweight[m_samples].second<<endl;

            int iRnd = mtrand.randInt(GoodPointSize-1);
            CXYVector<float> kV_Point(3);
            CXYMatrix<float> kM_MiddlePoints(m_iNumMiddleEndPoints,3);

            kV_Point = kMSamplesPoints.GetRow(GoodPointInd[iRnd]) + (*node);
            GetMiddlePoints((*node), kV_Point, kM_MiddlePoints);



            ptr[m_samples][j] -> append_child(node, kV_Point);



            for (int iNumPoint=m_iNumMiddleEndPoints-1; iNumPoint >= 0; iNumPoint --) {
                m_pOctree[m_samples] -> addPoint(kM_MiddlePoints[iNumPoint][0],kM_MiddlePoints[iNumPoint][1],kM_MiddlePoints[iNumPoint][2],1,m_maxDepth);
            }

            //cout<<"m_pOctree is "<<m_pOctree[j]<<endl;

        }
   //     ms_logweight[m_samples].first + = m_LogWeight[j].first;
       // sum = 0;

   //    cout<<"grow one chain is ended"<<endl;

        return Flag;

}

//----------------------------------------------------------------------GAMZE 10/9----------------------------------------------
void CXYEnsemble::Resampling(int k, int m_samples)
{

    int iNumNodes = GetNumNodes();

    vector <float> probability;

    float sumofweights = 0;

    vector <int> list;

    int  index=0;

    int listsize=0;

    int rand1;



    vector <vector <tree< CXYVector<float> >*> > ptr;

    cout<<" start fling"<<endl;

    int f=0;
   // int t = 0;
    // int t=k+(f*k);
    int m = iNumNodes;
    CXYVector<int> VectorofIndex;

    indicator=0;

    while (f==0 && GrowmChain(k,f,m_samples))
    {
     f=f+1;
    }

int t=2*k+(f*k);
       while (t<=iNumNodes && (GrowmChain(k,f,m_samples)))
        {
                double ESS = CalculateEffSampleSize(m_samples);
                cout<<"ESS is "<<ESS<<endl;
                if (ESS <= 0.3*m_samples)
                {

                    t=2*k+(f*k);
                    cout<<"t is "<<t<<endl;
                    VectorofIndex = getResampleIndex(m_samples);
                    CreateNewPointTree(t,m_samples, VectorofIndex);
                    cout<<"we created new point tree"<<endl;
                    CreateNewOctree(t,m_samples);
                    cout<<"we created new octree"<<endl;
                    ReconstructLogweight(m_samples, VectorofIndex);
                    cout<<"we reconstructed logweights"<<endl;
                   // if (t<=50)
                   // {
                        f=f+1;
                   // }
                   // else
                   // {
                   //     f=f+0.5;
                   // }
              }
              else
             {
                    f=f+1;
             }

            t=2*k+(f*k);


            }

    }

CXYVector <float> CXYEnsemble::getResampleprob(int m_sample)
{
    CXYVector<float> prob(m_sample);
    cout<< "prob for weights is being created"<<endl;
    vector <FloatIntPair> weights;

    for (int i = 0; i<m_sample; i++)
    {
        weights.push_back(FloatIntPair(0,1));
    }
    for (int i = 0; i<m_sample; i++)
    {
        weights[i].first = ms_logweight[i].first;
        weights[i].second = ms_logweight[i].second;
    }

    sort(weights.begin(), weights.end(), FloatIntPairCompare());

    float max_prob = weights[m_sample-1].first;
    cout<<"max prob is "<<max_prob<<endl;

    for (int i=0; i<m_sample; i++)
    {
        prob[i] = exp(ms_logweight[i].first-max_prob);
    }

    return prob;
}
CXYVector <int> CXYEnsemble::getResampleIndex(int m_sample)
{

    int m= m_sample;
   CXYVector<int> new_vector(m_sample);
    CXYVector<float> prob = getResampleprob(m_sample);
   // cout<<"are we here yet "<<endl;
    vector <FloatIntPair> p;
    for (int i = 0; i<m; i++)
    {
        p.push_back(FloatIntPair(0,1));
    }

    	
    for (int i = 0; i<m; i++)
    {
        p[i].first = prob[i];
        p[i].second =i;
    }

   	

    sort(p.begin(), p.end(), FloatIntPairCompare());
    	
    for (int i = 1; i<m; i++)
    {
        p[i].first = p[i].first+p[i-1].first;
      //  p[i].second =i;
    }


    float maxP= p[m-1].first;
   // cout<<"max weight is "<<maxWeight<<endl;
    cout<<"get index is working"<<endl;
    int s=0;
    //int j=0;

    while(s<m_sample)
    {
        float a_rnd=  static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/p[m-1].first));
      // cout<<"a_rnd is "<<a_rnd<<endl;
       // int j = rand() % m_sample;
      //  cout<<"j is "<<j<<endl;
      //  cout<<"prob[j] is "<<prob[j]<<endl;
      for (int j=0; j<m_sample; j++)
      {
        if (p[j].first > a_rnd && p[j].first>0)
        {
            new_vector[s]=p[j].second;
            s=s+1;
           // cout<<"j is "<<p[j].second<<endl;
            break;
        }
      }
    }
    cout<<"s is "<<s<<endl;
    return new_vector;


}
double CXYEnsemble::CalculateEffSampleSize(int m)
{
   // cout<<"are we here yet "<<endl;
    vector <FloatIntPair> weights;
    for (int i = 0; i<m; i++)
    {
        weights.push_back(FloatIntPair(0,1));
    }

    for (int i = 0; i<m; i++)
    {
        weights[i].first = ms_logweight[i].first;
        weights[i].second =ms_logweight[i].second;
    }

    sort(weights.begin(), weights.end(), FloatIntPairCompare());

    float maxWeight= weights[m-1].first;
    cout<<"max weight is "<<maxWeight<<endl;
   //  cout<<"max weight is "<<m_LogWeight[3].first<<endl;

    CXYVector <float> relativeWeight(m);
    double sumweight = 0;
    CXYVector <float> q(m);

    for (int i = 0; i < m; i++)
    {
      //  cout<<"are we here yet "<<endl;
    //    cout<<"weight is "<<m_LogWeight[i].first<<endl;
     //   cout<<"weight is "<<maxWeight<<endl;
        relativeWeight[i] = exp(ms_logweight[i].first-maxWeight);
        // relativeWeight[i] = 0;


        sumweight = sumweight + relativeWeight[i];

    }



    for (int i = 0; i < m; i++)
    {
        q[i] = relativeWeight[i] / sumweight;
    }

    double ESS = 0;

    for (int i = 0; i < m; i++)
    {
        ESS = ESS + pow(q[i],2);
    }

    ESS = 1 / ESS;

    return ESS;
}

void CXYEnsemble::CreateNewPointTree(int k, int m_samples, CXYVector<int> vectorofI)
{
    vector <vector < tree<CXYVector<float> >* > >  ptr = GetTree();
    vector < vector < CXYMatrix<float> > > newpoints;
    for (int i = 0; i < m_samples; i++)
    {

	     vector < CXYMatrix <float> >  v;
         for (int j=0; j<m_fNumofChains; j++)
         {

             CXYMatrix <float>  MatTemp;
             v.push_back(MatTemp);
            // delete TreeTemp;
         }
         newpoints.push_back(v);
        // v.clear();

	}



   // int iNumNodes= k;
    vector <  tree < CXYMatrix <float> >::iterator>  top, root;
    cout<<"starting the iteration in createnewpointtree"<<endl;
   for (int s=0; s<m_samples; s++)
   {

        for ( int j = 0; j < m_fNumofChains; j++ )
        {
            newpoints[s][j] = GetPoints(ptr[s][j],k) ;
          //  cout<<"get points is happening"<<endl;
        }

   }

    cout<<"get points ended"<<endl;





   for (int s=0; s<m_samples; s++)
   {

        vector <  tree < CXYVector <float> >::iterator>  top, root;
     //   cout<<"index is "<<ms_logweight[m_samples-s-1].second<<endl;
        CXYMatrix<float> newnewpoints;

        for (int i=0; i<m_fNumofChains; i++)
        {
            newnewpoints = newpoints[vectorofI[s]][i];

    //     if ( !ptr.at(m_LogWeight[j].second) -> empty() )
    //        {
            ptr[s][i] -> clear();
     //       }
            top.push_back(ptr[s][i] -> begin());
            float fCoord[3] = {newnewpoints[0][0],newnewpoints[0][1],newnewpoints[0][2]};
            CXYVector<float> kV_startpoint (3,fCoord);
            root.push_back((ptr[s][i] -> insert (top.at(i), kV_startpoint)));


            for ( int j = 1; j < k; j++ )
            {

                tree< CXYVector<float> >::pre_order_iterator node;

                node = ptr[s][i] -> end();
                node --;
                CXYVector<float> kV_Point(3);

            //cout<<"i is "<<i<<endl;

                 kV_Point[0] = newnewpoints[j][0];
                 kV_Point[1] = newnewpoints[j][1];
                 kV_Point[2] = newnewpoints[j][2];


                 ptr[s][i] -> append_child(node, kV_Point);
           //  cout<<"is this happening"<<endl;
           //cout<<"i is "<<i<<endl;
         //  cout<<"iNumnodes is "<<k<<endl;
            }
        }

      //  newnewpoints.Deallocate();
   }


cout<<"create new points ended"<<endl;

    //newpoints.clear();

   // newpoints.clear();
}

void CXYEnsemble::CreateNewOctree(int k, int m_samples)

{


    vector <vector <tree<CXYVector<float> >* > > ptr = GetTree();

    OcTree<float,int>* p_octree;

      for ( int i = 0; i < m_pOctree.size(); i++ )
        {

            //if (m_pOctree[i] -> root_) {
               // cout<<"is this hapening"<<endl;
                delete m_pOctree[i];

            //}
        }


    m_pOctree.clear();

    //int iNumNodes= k;
        for (int i = 0; i < m_samples; i++)
	{
        boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  =  boost::make_shared< OcTreeNodeAllocator< float , int > >();
        m_pOctree.push_back(new OcTree<float,int>(m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator) );

	}


    for( int j=0; j<m_samples; j++)
    {

 //       a = GetPoints(ptr[j], iNumNodes, j);

        for (int s=0; s<m_fNumofChains; s++)
        {

        CXYMatrix<float> new_points(k,3) ;
        new_points = GetPoints(ptr[j][s], k);



         //   cout<<"new_points are "<<new_points[5][1]<<endl;
     //       boost::shared_ptr< OcTreeNodeAllocator< float, int > > allocator  =  boost::make_shared< OcTreeNodeAllocator< float , int > >();
    // p_octree = new OcTree<float,int> (m_center,m_minimumVolumeSize,m_fNucleusSphereDiameter,allocator);
    // m_pOctree.push_back(p_octree);

     m_pOctree[j] ->addPoint(new_points[0][1],new_points[0][2],new_points[0][3],1,m_maxDepth);

        for (int i=1; i<k; i++)

        {


            CXYMatrix < float > kMSamplesPoints = GetSamplesOrg();
            kMSamplesPoints *= GetSegLength(i-1);
            float fCoord_pre[3];
            float fCoord[3];
            CXYVector<float> kV_Point(3);
            CXYVector<float> kV_Point_pre(3);
            CXYMatrix<float> kM_MiddlePoints(m_iNumMiddleEndPoints,3);

            for (int index=0; index<3; index++)
            {
                fCoord_pre[index]=new_points[i-1][index];
                fCoord[index] = new_points[i][index];
               // cout<<"kv points are "<<fCoord[index]<<endl;
            }


            kV_Point = CXYVector<float> (3,fCoord);
            kV_Point_pre = CXYVector<float> (3,fCoord_pre);

            GetMiddlePoints(kV_Point_pre, kV_Point, kM_MiddlePoints);

            for (int iNumPoint=m_iNumMiddleEndPoints-1; iNumPoint >= 0; iNumPoint --) {
                m_pOctree[j] -> addPoint(kM_MiddlePoints[iNumPoint][0],kM_MiddlePoints[iNumPoint][1],kM_MiddlePoints[iNumPoint][2],1,m_maxDepth);
            }

                    //    cout<<"is this happening"<<endl;
              //      new_points.Deallocate();
                          //  delete p_octree;
        }

        }

        //cout<<"is this happening"<<endl;
   //     cout<<"j is "<<j<<endl;
       // newpoints.clear();
        //delete[] newnewpoints;



    }


}

CXYMatrix<float>  CXYEnsemble::GetPoints(tree< CXYVector<float> >*  ptr, int k)
{

    tree < CXYVector <float> > ::iterator pos;
    CXYMatrix<float> oldpoints(k,3);
 //   vector < CXYMatrix<float> > a;

 //   for (int s=0; s<m; s++)
  //  {
    pos = ptr -> begin();
    int i=0;
    while (ptr-> is_valid(pos))
    {
     //   cout<<"is this happening"<<endl;
        oldpoints[i][0]=(*pos)[0];
        oldpoints[i][1]=(*pos)[1];
        oldpoints[i][2]=(*pos)[2];
    //    cout<<"oldpoints are "<<oldpoints[i][2]<<endl;
        ++pos;
        i=i+1;
      //  cout<<"i is "<<i<<endl;
    }

   return oldpoints;

 //   }

 //  cout<<"get points is ended"<<endl;

}

void CXYEnsemble::ReconstructLogweight(int m, CXYVector<int> vectorofI)
{
   CXYVector<double> newweights(m);
    CXYVector<double> temp_weight(m);
    CXYVector<float> prob = getResampleprob(m);
    float sumofweights=0;
    float sumoftemp = 0;

    for (int i=0; i<m; i++)
    {
        newweights[i] = ms_logweight[i].first;
        sumofweights=sumofweights+newweights[i];
    }




    for (int i=0; i<m; i++)
    {
        temp_weight[i] = newweights[vectorofI[i]]-(log(prob[vectorofI[i]]));
        sumoftemp=sumoftemp+temp_weight[i];
      //  cout<<"log of prob is "<<log(2*prob[VectorofIndex[i]])<<endl;
    }

    float alpha = (sumofweights-sumoftemp)/m;

    for (int i=0; i<m; i++)
    {
	ms_logweight[i].first = alpha + temp_weight[i];
    }

   // newweights.Deallocate();
}
//---------------------here I am---------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void CXYEnsemble::WritemChain(char* fn, int m, int m_samples)
{
//http://stdcxx.apache.org/doc/stdlibug/34-2.html

	vector <vector < tree<CXYVector<float> >* > >  ptr = GetTree();




// int count =0;


	char buffer[1024];

for (int j=0; j<m_samples; j++)
{

    ostream* fp;
        if (strcmp(fn,"")) {
        sprintf(buffer, "%s/%s_%d", GetOutPath(),fn,j);
        cout<<buffer<<endl;
		fp = new std::ofstream(buffer);
        }else {
            fp = &std::cout;
        }



       *fp << "# LogWeight= " << ms_logweight[j].first << endl;


	for (int i=0; i < m_fNumofChains; i++)

	{
	    int count = 0;
	     tree < CXYVector <float> >::iterator pos;

	    *fp<< "#  Chr "<<(i+1)<<endl;

	    pos=ptr[j][i]->begin();
        //cout<<"chain lengths are "<<m_fChainLengths[i]<<endl;
    //    while (ptr[j][i]->is_valid(pos)) {
         while (count<m_fChainLengths[i]){

            sprintf(buffer,"%3e\t%3e\t%3e",(*pos)[0],(*pos)[1],(*pos)[2]);
        //cout<<m_LogWeight[m_LogWeight[i].second].first<<endl;


            *fp << buffer <<endl;
		//cout<<buffer<<endl;
            ++pos;
            count=count+1;
		//cout<<"pos is "<<count<<endl;
        }

	}

	if (fp != &std::cout)
		delete fp;

	}



    cout<< "is this happening" <<endl;

}
//----------------------------------------------------------------------------

void CXYEnsemble::SetOutPath(char* cPathName)
{
	CXYFile::MakeDirectory(cPathName,0755);
  cout << "make directory " << cPathName << endl;
	strcpy(m_cOutPath,cPathName);
}

//----------------------------------------------------------------------------
char* CXYEnsemble::GetOutPath(void)
{
	return m_cOutPath;
}
//----------------------------------------------------------------------------
// void CXYEnsemble::SetSegLengths(char* cStartEndFile){
// 	CXYMatrix<float> kMStartEnd = CXYFile::ReadMatrix(cStartEndFile);
// 	for (int i= 0; i<kMStartEnd.GetRows(); i++) {
// 		int baselen = (int) ( kMStartEnd[i][1] - kMStartEnd[i][0] );
// 		m_vfSegLength.push_back(((float)baselen * GetPackingDensity()));
// 	}
// }
//----------------------------------------------------------------------------
void CXYEnsemble::SetSegLengths(char* cStartEndFile,const char* cMethod){
	CXYMatrix<float> kMStartEnd = CXYFile::ReadMatrix(cStartEndFile);

	if (strcmp(cMethod, "AVG") == 0) {
		// average different node length
		float len = 0L;
		for (int i= 0; i<kMStartEnd.GetRows(); i++) {
			len = len + ( kMStartEnd[i][1] - kMStartEnd[i][0] + 1 ) * GetPackingDensity();
		}
		float avglen = len/(float)kMStartEnd.GetRows();
		// set fixed length
		SetPersistenceLength( avglen );

		for (int i= 0; i<kMStartEnd.GetRows(); i++) {
			m_vfSegLength.push_back(GetPersistenceLength());
		}

	} else if (strcmp(cMethod, "DIF") == 0) {
		// different node length
		for (int i= 0; i<kMStartEnd.GetRows(); i++) {
			int baselen = (int) ( kMStartEnd[i][1] - kMStartEnd[i][0] );
			m_vfSegLength.push_back(((float)baselen * GetPackingDensity()));
		}
	} else {
		cerr	<< "Please type Method" << endl;
		exit(-1);
	}
}
//----------------------------------------------------------------------------

vector<float> & CXYEnsemble::GetSegLengths(void){
	return m_vfSegLength;
}
//----------------------------------------------------------------------------
float CXYEnsemble::GetSegLength(int ind){
	return m_vfSegLength[ind];
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetSegLengths(void){
	for (int i= 0; i<GetNumNodes(); i++) {
		m_vfSegLength.push_back(GetPersistenceLength());
	}
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetContIndex(char* cContIndFile){
//	CXYVector<int> vtmp = CXYFile::ReadVectorInt(cContIndFile);
//	for (int i =0; i<vtmp.GetSize(); i++) {
//		// index start from 0
//		m_MContInd.push_back(vtmp[i]-1);
//	}
  m_MContInd = CXYFile::ReadMatrix(cContIndFile);
}
//----------------------------------------------------------------------------
void CXYEnsemble::SetContIndex(void){
  m_MContInd.SetSize(m_iNumNodes,3);
  for (int i=0; i<m_iNumNodes; i++) {
    m_MContInd[i][0] = i+1;
    m_MContInd[i][1] = i+1;
    m_MContInd[i][2] = i+1;
  }
}
//----------------------------------------------------------------------------
void CXYEnsemble::GetMiddlePoints(CXYVector<float>& StartPoint, CXYVector<float>& EndPoint, CXYMatrix<float> & MiddleEndPoints){

  float dx = EndPoint[0] - StartPoint[0];
  float dy = EndPoint[1] - StartPoint[1];
  float dz = EndPoint[2] - StartPoint[2];

  for (int i =0; i< m_iNumMiddleEndPoints; i++) {
    MiddleEndPoints[i][0] = StartPoint[0] + ((i+1)/(float)m_iNumMiddleEndPoints) * dx;
    MiddleEndPoints[i][1] = StartPoint[1] + ((i+1)/(float)m_iNumMiddleEndPoints) * dy;
    MiddleEndPoints[i][2] = StartPoint[2] + ((i+1)/(float)m_iNumMiddleEndPoints) * dz;
  }

}
//----------------------------------------------------------------------------
bool CXYEnsemble::IsSatisfyCondition(CXYMatrix<float> & MiddleEndPoints, int j){
  bool flag = true;
  for (int i=0; i< m_iNumMiddleEndPoints; i++) {
     //    cout<< "nasi lan "<<MiddleEndPoints.GetRow(i)<<endl;



    if (! IsInsideSphere(MiddleEndPoints.GetRow(i))  ) {
        //cout << "i is "<<i<<endl;

      flag = false;
      break;
    }

    //cout<<"WTF1"<<endl;

    if (IsCollision(MiddleEndPoints.GetRow(i),j)){
//        cout<<"this has happenned"<<endl;
 //     cout << MiddleEndPoints[m_iNumMiddleEndPoints-1][0] << " " << MiddleEndPoints[m_iNumMiddleEndPoints-1][1] << " " << MiddleEndPoints[m_iNumMiddleEndPoints-1][2] << endl;
      flag = false;
      break;
    }
  }
  //cout<<"flag is"<<flag<<endl;
  return flag;
}
//----------------------------------------------------------------------------
bool CXYEnsemble::IsCollision(CXYMatrix<float> & MiddleEndPoints, int j){
 //   cout<<"WTF"<<endl;
  bool flag = false;

  for (int i=0; i< m_iNumMiddleEndPoints; i++) {
  //    cout<<"this is is "<<i<<endl;
    if ( IsCollision(MiddleEndPoints.GetRow(i),j)){
          //  cout << MiddleEndPoints[m_iNumMiddleEndPoints-1][0] << " " << MiddleEndPoints[m_iNumMiddleEndPoints-1][1] << " " << MiddleEndPoints[m_iNumMiddleEndPoints-1][2] << endl;
      flag = true;
      break;
    }
  }
  return flag;
}

//----------------------------------------------------------------------------
void CXYEnsemble::GetGoodPoints( CXYVector<float>& prvnode, CXYMatrix<float>& kMSamplesPoints,vector<int>& GoodPointInd, vector<int>& NoCollisionPointInd, int j){

  //cout << prvnode[0] << " " << prvnode[1] << " " << prvnode[2] << endl;

  CXYVector<float> kV_Point(3);
  CXYMatrix<float> kM_MiddlePoints(m_iNumMiddleEndPoints,3);

  int iSampleSize = kMSamplesPoints.GetRows();


//  cout << "===================="<< iSampleSize << endl;
  for (int i = 0; i < iSampleSize; i++) {
    kV_Point = kMSamplesPoints.GetRow(i) + prvnode;

   // cout << kV_Point[0] << " " << kV_Point[1] << " " << kV_Point[2] << endl;

    GetMiddlePoints(prvnode, kV_Point, kM_MiddlePoints);

 //   cout << "i = " << i <<endl;
 //   cout<<"km middle points is "<<kM_MiddlePoints[1][2]<<endl;
 //   cout<<"j is "<<j<<endl;
    if (IsSatisfyCondition(kM_MiddlePoints, j)){
    //    cout << "this part has happenned" <<endl;
      GoodPointInd.push_back(i);
    }
    // No collision
    if (! IsCollision(kM_MiddlePoints, j)){
      NoCollisionPointInd.push_back(i);
    }
  }
//  cout << "===================="<< iSampleSize << endl;
    //cout<<"is this happening"<<endl;


}

//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
