#include <stdio.h>
#include <iostream>
#include <fstream>
#include <list>
#include <math.h>
#include <cstring>
#include <time.h>
#include <sstream>
#include "NOrMAL.h"

bool compare_int (int i,int j) { return (i<j); }

Nucs::Nucs(){
	_sigma = 15;
	_delta_f = 0;
	_delta_r = 0;
	_nuc_area = 120; //area occupied by one nucleosome and it's linker
	_nuc_size_prior =  73; // 1/2 size of the nucleosome size
	_prior_var = 0.001;
    _prior_end = 0.01;
	_eps = 0;
	_shift_max =  100;
	_shift_min = 35;
	_sigma_max = 60;
	_sigma_min = 3;
	_merge_distance = 0.15;
	_read_size = 0 ;
}

double Nucs::PriorEnd() 
{
    return _prior_end;
}

int Nucs::NucSizePrior()
{
    return _nuc_size_prior;
}

void Nucs::PriorSet(double x)
{
    _prior_var = x;
}

int Nucs::Refine()
{
    int flag = 0; // counter of how many mergings were performed
    std::list<Nuc_Model>::iterator it, itl, itr;
    it = Param.begin();
    it++;
    while( it != Param.end() )
    {
        itl = it;
        itl--;
        itr = it;
        itr++;
        int distanceL = -it->mju+it->shift + itl->mju + itl->shift;
        int distanceR = it->mju + it->shift - itr->mju + itr->shift;
        if ( ( distanceL > _merge_distance*itl->shift ) || ( distanceR > _merge_distance*itr->shift ) )
        {
            flag++;
            //merging two nucleosomes
            if( distanceL <= distanceR )
            {
                double prob_new = it->prob + itl->prob;
                int mju_new = ( itl->mju * itl->prob + it->mju * it->prob )/prob_new;
                double shift_new =  ( itl->shift * itl->prob + it->shift * it->prob )/prob_new;
                double sigma_new = sqrt( ( pow( itl->sigma, 2) * itl->prob + pow( it->sigma, 2) * it->prob )/prob_new );
                double delta_f_new = 0;//sqrt( ( pow( itl->delta_f, 2) * itl->prob + pow( it->delta_f, 2) * it->prob )/prob_new );
                double delta_r_new = 0;//sqrt( ( pow( itl->delta_r, 2) * itl->prob + pow( it->delta_r, 2) * it->prob )/prob_new );
                itl->mju = mju_new;
                itl->sigma = sigma_new;
                itl->shift = shift_new;
                itl->delta_f = delta_f_new;
                itl->delta_r = delta_r_new;
                itl->prob = prob_new;
                it = Param.erase(it);
            }
            else
            {
                double prob_new = it->prob + itr->prob;
                int mju_new = ( itr->mju * itr->prob + it->mju * it->prob )/prob_new;
                double shift_new =  ( itr->shift * itr->prob + it->shift * it->prob )/prob_new;
                double sigma_new = sqrt( ( pow( itr->sigma, 2) * itr->prob + pow( it->sigma, 2) * it->prob )/prob_new );
                double delta_f_new = 0;//sqrt( ( pow( itl->delta_f, 2) * itl->prob + pow( it->delta_f, 2) * it->prob )/prob_new );
                double delta_r_new = 0;//sqrt( ( pow( itl->delta_r, 2) * itl->prob + pow( it->delta_r, 2) * it->prob )/prob_new );
                it->mju = mju_new;
                it->sigma = sigma_new;
                it->shift = shift_new;
                it->delta_f = delta_f_new;
                it->delta_r = delta_r_new;
                it->prob = prob_new;
                it = Param.erase(itr);
            }
        }
        else
        {
            it++;
        }
    }
    return flag;
}

int Nucs::Load_Data_X(char* str) //Load datapoint from txt file, forward reads
{
    std::ifstream f_in;
    f_in.open( str );
    int x;
    int id_count = 0;
    if( !f_in.is_open() )
    {
        std::cerr<<"Can't open input file for forward reads"<<std::endl;
        return 0;
    }

    while( !f_in.eof() )
    {
        id_count++;
        f_in >> x;
        X.push_back( x );
    }
    f_in.close();
    X.sort();   
    X.pop_back(); 
    return 1;
}

int Nucs::Load_Data_Y(char* str) //Load datapoint from txt file, reverse reads
{
    std::ifstream f_in;
    f_in.open( str );
    int x;
    int id_count = 0;
    if( !f_in.is_open() )
    {
        std::cerr<<"Can't open input file for reverse reads"<<std::endl;
        return 0;
    }
    while( !f_in.eof() )
    {
        id_count++;
        f_in >> x;
        Y.push_back( x+_read_size );
    }
    f_in.close();
    Y.sort();   
    Y.pop_back(); 
    return 1;
}
     
int Nucs::Init(int start, int finish, int K) // Initialize nucleosome positions
{
    int region_size = finish - start;
    double linker = region_size / (K+1.0);
    Nuc_Model P;
    for (int i=0; i<=K+3; i++)
    {
        P.mju = start + (i-1) * linker ;
        P.sigma = _sigma;
        P.delta_f = _delta_f;
        P.delta_r = _delta_r;
        P.shift = _nuc_size_prior;
        P.prob = 1.0/K;
        P.prob_f = 0;
        P.prob_r = 0;
        Param.push_back( P );
    }
    Param.front().prob = 0;
    Param.back().prob = 0;
    return 1;
}

int Nucs::InitTF( char* file )
{
    std::ifstream input;
    input.open(file);
    if( input.is_open() )
    {
        std::string line;
        getline(input,line);
        while( getline(input,line) )
        {
            
            std::stringstream str;
            str<< line;
            std::string chr;
            int center, width,  numForwardReads, numReverseReads;
            double corr_score;
            str>> chr >> center >> width >> corr_score >> numForwardReads >> numReverseReads;  
            Nuc_Model P;
            P.mju = center ;
            P.sigma = _sigma;
            P.delta_f = _delta_f;
            P.delta_r = _delta_r;
            P.shift = width/2;
            P.prob = 1;
            P.prob_f = 0;
            P.prob_r = 0;
            Param.push_front( P );
        }
        input.close();
        return 1;
    }
    else
    {
        std::cerr<<"No initialization files"<<std::endl;
        return 0;
    }
}

int Nucs::Print( char* str)
{
    std::ofstream f_out;
    f_out.open( str );
    std::list<Nuc_Model>::iterator it, B, E;
    B = Param.begin();
    E = Param.end();
    B++;
    E--;
    f_out<< "<Position> <Fuzziness> <Size> <Confidence Score> <Forward votes> <Reverse votes> "  <<  std::endl;
    for ( it = B; it != E; it++ )
    {
        if( it->shift > 0 && it->shift < 250) //filtering abnormal size clusters 
            f_out<< it->mju << " " << (int)it->sigma <<" "<<2*it->shift << " " << it->prob << " " << it->prob_f << " " << it->prob_r <<  std::endl;
    }
    f_out.close();    
    return 1;
}

inline double F_pdf(int x, double mju, double sigma2)
{
    return exp( - pow( x - mju + .0, 2) / sigma2 ) / sqrt(2*M_PI*sigma2);
}


int Nucs::hardEM_update()
{
    int counter_delete = 0; //counter of how many nucleosomes were deleted
    std::list<Nuc_Model> P;
    P = Param;
    std::list<Temp_Param> counters;
    Temp_Param T;
    for( unsigned int i=0; i<P.size();i++)
    {
        T.X_hat = 0;
        T.Y_hat = 0;
        T.X2_hat = 0;
        T.Y2_hat = 0;
        T.Tij_X = 0;
        T.Tij_Y = 0;
        T.Dx_hat = 0;
        T.Dy_hat = 0;
        counters.push_back(T);
    }
    int i_d=0;
    std::list<int>::iterator it;
    std::list<Nuc_Model>::iterator jt, jt0, jtl, jtr, b;
    std::list<Temp_Param>::iterator cur,cur0;
    b = P.end();
    b--;
    jt0 = P.begin();
    jt0++;
    cur0 = counters.begin();
    cur0++;
    bool flag_processed;
    i_d = 0;
    for (it = X.begin(); it!=X.end(); it++)
    {   
        flag_processed = false;
        int Xi = *it;
        for( jt = jt0, cur = cur0; jt!=b; jt++, cur++)
        {
            i_d++;
            jtl = jt;
            jtl--;
            jtr = jt;
            jtr++;
            int distance = fabs(Xi - jt->mju + jt->shift);
            if( ( distance <= fabs(Xi - jtl->mju + jtl->shift) ) && ( distance <= fabs(Xi - jtr->mju + jtr->shift) ) )
            {
                flag_processed = true;
                jt0 = jtl;
                cur0 = cur;
                cur0--;
                //all computations go here
                double temp=0;
                double Tij = 0;
                temp = jt->prob*F_pdf( Xi, jt->mju - jt->shift, pow(jt->sigma, 2)+pow(jt->delta_f, 2) );
                Tij=temp;
                    cur->X_hat += Tij * Xi;
                    cur->X2_hat += Tij * pow( Xi + jt->shift - jt->mju, 2);
                    cur->Dx_hat += Tij * jt->shift;
                    cur->Tij_X += Tij;
            }
            else
            {
                if( flag_processed )
                {
                    i_d--;
                    break;
                }
            }
        }
    }

    //do the same for Y (reverse tags)
    jt0 = P.begin();
    jt0++;
    cur0 = counters.begin();
    cur0++;
    i_d = 0;
    for (it = Y.begin(); it!=Y.end(); it++)
    {   
        flag_processed = false;
        int Xi = *it;
        for( jt = jt0, cur = cur0; jt!=b; jt++, cur++)
        {
            i_d++;
            jtl = jt;
            jtl--;
            jtr = jt;
            jtr++;
            int distance = fabs(Xi - jt->mju - jt->shift);
            if( ( distance <= (Xi - jtl->mju - jtl->shift) ) && ( distance <= fabs(Xi - jtr->mju - jtr->shift) ) )
            {
                flag_processed = true;
                jt0 = jtl;
                cur0 = cur;
                cur0--;
                //all computations go here
                double temp=0;
                double Tij = 0;
                temp = jt->prob*F_pdf( Xi, jt->mju + jt->shift, pow(jt->sigma, 2)+pow(jt->delta_f, 2) );
                Tij=temp;
                    cur->Y_hat += Tij * Xi;
                    cur->Y2_hat += Tij * pow( Xi - jt->shift - jt->mju, 2);
                    cur->Dy_hat += Tij * jt->shift;
                    cur->Tij_Y += Tij;
            }
            else
            {
                if( flag_processed )
                {
                    i_d--;
                    break;
                }
            }
        }
    }
    //...//


    double norm_const = 0;
    jt0 = Param.begin();
    jt0++;
    cur0 = counters.begin();
    cur0++;
    b = Param.end();
    b--;
    
    //update rules go here...
    jt = jt0;
    cur = cur0;
    i_d = 0;
    while ( jt != b )
    {
        i_d++;
        jt->prob_f = cur->Tij_X;
        jt->prob_r = cur->Tij_Y;
        double temp;
        temp = cur->Tij_X + cur->Tij_Y;
        if ( temp > _eps )
        {
            jt->mju = (int) ( ( cur->X_hat + cur->Dx_hat + cur->Y_hat - cur->Dy_hat ) / temp ) ;
            if ( cur->Tij_X > 0 )
            {
                jt->delta_f = sqrt( cur->X2_hat / cur->Tij_X);// - pow(jt->sigma,2) );
            }
            else
            {
                jt->delta_f = 0;
            }
            if ( cur->Tij_Y > 0 )
            {
                jt->delta_r = sqrt( cur->Y2_hat / cur->Tij_Y);// - pow(jt->sigma,2) );
            }
            else
            {
                jt->delta_r = 0;
            }
            jt->sigma = sqrt( ( cur->X2_hat + cur->Y2_hat ) / temp );

            if(( cur->Tij_X > 0 )&&( cur->Tij_Y > 0 ))
            {
                jt->shift = ( -cur->X_hat/cur->Tij_X + cur->Y_hat/cur->Tij_Y  )/2 ;
            }
            else
                jt->shift = 0; 
            norm_const += jt->prob; 
            cur++;
            jt++;
        }
        else
        {
            counter_delete++;
            jt = Param.erase(jt);
            cur++;
        }
    }
    int mult = Param.size();
    jt = Param.begin();
    while( jt!= b )
    {
        jt->prob = jt->prob * mult / norm_const;
        jt++;
    }

    return counter_delete;
}

int Nucs::EM_update()
{
    int counter_delete = 0; //counter of how many nucleosomes were deleted
    std::list<Nuc_Model> P;
    P = Param;
    std::list<Temp_Param> counters;
    Temp_Param T;
    for( unsigned int i=0; i<P.size();i++)
    {
        T.X_hat = 0;
        T.Y_hat = 0;
        T.X2_hat = 0;
        T.Y2_hat = 0;
        T.Tij_X = 0;
        T.Tij_Y = 0;
        T.Dx_hat = 0;
        T.Dy_hat = 0;
        counters.push_back(T);
    }
    int i_d=0;
    std::list<int>::iterator it;
    std::list<Nuc_Model>::iterator jt, jt0, jtl, jtr, b;
    std::list<Temp_Param>::iterator cur,cur0;
    b = P.end();
    b--;
    jt0 = P.begin();
    jt0++;
    cur0 = counters.begin();
    cur0++;
    bool flag_processed;
    i_d = 0;
    for (it = X.begin(); it!=X.end(); it++)
    {   
        flag_processed = false;
        int Xi = *it;
        for( jt = jt0, cur = cur0; jt!=b; jt++, cur++)
        {
            i_d++;
            jtl = jt;
            jtl--;
            jtr = jt;
            jtr++;
            if( ( Xi >= jtl->mju ) && ( Xi <= jtr->mju ) )
            {
                flag_processed = true;
                jt0 = jtl;
                cur0 = cur;
                cur0--;
                //all computations go here
                double temp=0;
                double Tij = 0;
                temp = jt->prob*F_pdf( Xi, jt->mju - jt->shift, pow(jt->sigma, 2)+pow(jt->delta_f, 2) );
                Tij = temp / ( temp + jtl->prob*F_pdf( Xi, jtl->mju - jtl->shift, pow(jtl->sigma, 2)+pow(jtl->delta_f, 2) ) + jtr->prob*F_pdf( Xi, jtr->mju - jtr->shift, pow(jtr->sigma, 2)+pow(jtr->delta_f, 2) ) );
                if( isnan( Tij ) ) 
                {
                    Tij = 0;
                }
                else
                {
                    cur->X_hat += Tij * Xi;
                    cur->X2_hat += Tij * pow( Xi + jt->shift - jt->mju, 2);
                    cur->Dx_hat += Tij * jt->shift;
                    cur->Tij_X += Tij;
                }
            }
            else
            {
                if( flag_processed )
                {
                    i_d--;
                    break;
                }
            }
        }
    }

    //do the same for Y
    jt0 = P.begin();
    jt0++;
    cur0 = counters.begin();
    cur0++;
    i_d = 0;
    for (it = Y.begin(); it!=Y.end(); it++)
    {   
        flag_processed = false;
        int Xi = *it;
        for( jt = jt0, cur = cur0; jt!=b; jt++, cur++)
        {
            i_d++;
            jtl = jt;
            jtl--;
            jtr = jt;
            jtr++;
            if( ( Xi >= jtl->mju ) && ( Xi <= jtr->mju ) )
            {
                flag_processed = true;
                jt0 = jtl;
                cur0 = cur;
                cur0--;
                //all computations go here
                double temp=0;
                double Tij = 0;
                temp = jt->prob*F_pdf( Xi, jt->mju + jt->shift, pow(jt->sigma, 2)+pow(jt->delta_r, 2) );
                Tij = temp / ( temp + jtl->prob*F_pdf( Xi, jtl->mju + jtl->shift, pow(jtl->sigma, 2)+pow(jtl->delta_r, 2) ) + jtr->prob*F_pdf( Xi, jtr->mju + jtr->shift, pow(jtr->sigma, 2)+pow(jtr->delta_r, 2) ) );
                if ( isnan( Tij ) ) 
                {
                    Tij = 0;
                }
                else
                {
                    cur->Y_hat += Tij * Xi;
                    cur->Y2_hat += Tij * pow( Xi - jt->shift - jt->mju, 2);
                    cur->Dy_hat += Tij * jt->shift;
                    cur->Tij_Y += Tij;
                }
            }
            else
            {
                if( flag_processed )
                {
                    i_d--;
                    break;
                }
            }
        }
    }
    //...//


    double norm_const = 0;
    jt0 = Param.begin();
    jt0++;
    cur0 = counters.begin();
    cur0++;
    b = Param.end();
    b--;
    //update rules go here...
    jt = jt0;
    cur = cur0;
    i_d = 0;
    int _X_size = X.size();
    int _Y_size = Y.size();
    while ( jt != b )
    {
        i_d++;
        double _Tij_X = cur->Tij_X;
        double _Tij_Y = cur->Tij_Y;
        double _X_hat = cur->X_hat;
        double _Y_hat = cur->Y_hat;
        double _Dx_hat = cur->Dx_hat;  
        double _Dy_hat = cur->Dy_hat; 
        double _X2_hat = cur->X2_hat; 
        double _Y2_hat = cur->Y2_hat; 
        jt->prob_f = _Tij_X; 
        jt->prob_r = _Tij_Y;
        double temp;
        temp = _Tij_X + _Tij_Y;
        if ( temp > _eps )
        {
            jt->mju = (int) ( ( _X_hat + _Dx_hat + _Y_hat - _Dy_hat ) / temp ) ;
//           if ( cur->Tij_X > 0 )
//            {
//                jt->delta_f = sqrt( cur->X2_hat / cur->Tij_X - pow(jt->sigma,2) );
//            }
//            else
//            {
                jt->delta_f = 0;
//            }
//            if ( cur->Tij_Y > 0 )
//            {
//                jt->delta_r = sqrt( cur->Y2_hat / cur->Tij_Y - pow(jt->sigma,2) );
//            }
//            else
//            {
                jt->delta_r = 0;
//            }

//            if ( jt->delta_f > _sigma_max)
//                jt->delta_f = _sigma_max;
//            if ( jt->delta_f < _sigma_min )
//                jt->delta_f = _sigma_min;

//            if ( jt->delta_r > _sigma_max)
//                jt->delta_r = _sigma_max;
//            if ( jt->delta_r < _sigma_min )
//                jt->delta_r = _sigma_min;
            jt->sigma = sqrt( ( _X2_hat + _Y2_hat ) / temp );

            if ( jt->sigma > _sigma_max)
                jt->sigma = _sigma_max;
            if ( jt->sigma < _sigma_min )
                jt->sigma = _sigma_min;
            double _jt_sigma = jt->sigma; 
            jt->shift = ( -_X_hat + _Y_hat + _nuc_size_prior * pow(_jt_sigma / _prior_var, 2) ) / ( temp + pow(_jt_sigma / _prior_var, 2) ) ;
            double _jt_shift = jt->shift;
            if ( _jt_shift > _shift_max)
                jt->shift = _shift_max;
            if ( _jt_shift < _shift_min)
                jt->shift = _shift_min;

            jt->prob = temp / ( _X_size + _Y_size);
            norm_const += jt->prob; 
            cur++;
            jt++;
        }
        else
        {
            counter_delete++;
            jt = Param.erase(jt);
            cur++;
        }
    }

    int mult = Param.size();
    jt = Param.begin();
    while( jt!= b )
    {
        jt->prob = jt->prob * mult / norm_const;
        jt++;
    }

    return counter_delete;
}

int Nucs::Config(char* config)
{
    std::ifstream conf (config, std::ifstream::in);
    char line[256];
    std::string var;
    double val;
    while( conf.good() )
    {
        conf.getline( line, 256);
        if ( line[0] != '#' )
        {
            std::stringstream temp(line);
            temp >> var;
            temp >> val;
            if (var == "m_d")
                _merge_distance = val;
            else if (var == "n_s")
                _nuc_size_prior = val /2.0;
            else if (var == "n_a")
                _nuc_area = val;
            else if (var == "s_max")
                _shift_max = val /2.0;
            else if (var == "s_min")
                _shift_min = val /2.0;
            else if (var == "v_max")
                _sigma_max = val;
            else if (var == "read_size")
                _read_size = val;
        }
    }
    return 1;
}
