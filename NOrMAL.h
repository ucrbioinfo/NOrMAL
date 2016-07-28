#include <list>

struct Nuc_Model // Description of parameters to model nucleosomes
{
    public:
        int mju; // position of a nucleosome
        double sigma; // nucleosome fuzziness 
        double prob; // probability/score of a nucleosome
        double delta_f;
        double delta_r;
        double shift; // 1/2*size of the nucleosome
        double prob_f;
        double prob_r;
        int id; //unique ID of a nucleosome
};

struct Temp_Param //temperary parameter to store statistics
{
    public:
        double X_hat;
        double Y_hat;
        double X2_hat;
        double Y2_hat;
        double Tij_X;
        double Tij_Y;
        double Dx_hat;
        double Dy_hat;
};

class Nucs{
	public:
        Nucs();
		std::list<int> X,Y; //list of "forward" and "reverse" points respectevly
		std::list<Nuc_Model> Param; // list of model parameters, each node is separate cluster
		int Refine();
		int Load_Data_X(char* str); //Load datapoint from txt file to X list
		int Load_Data_Y(char* str); //Load datapoint from txt file to Y list
		int Init( int start, int finish, int K); // Initialize k nucleosome positions equally spaced in [start,finish] range
        int InitTF( char* file ); //Initialize nucleosomes using "file" as starting point, the format is output of Template Filtering tool
		int Print( char* str); //Print current set of Model parameters to file "str"
		int EM_update(); //Expectation Maximization (EM) iteration
        int hardEM_update(); //hard EM iteration
        int Config(char* config);
        void PriorSet(double x);
        double PriorEnd(); 
        int NucSizePrior();
		int _read_size ;
    private:
		double _sigma ;
		double _delta_f;
		double _delta_r ;
		double _nuc_area ; //area occupied by one nucleosome and it's linker
		double _nuc_size_prior ; // size of the nucleosome size
		double _prior_var ;
        double _prior_end;
		double _eps;
		double _shift_max;
		double _shift_min;
		double _sigma_max;
		double _sigma_min;
		double _merge_distance;
};
