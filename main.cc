#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <time.h>
#include <sstream>
#include "NOrMAL.h"

using namespace std;

int main(int argc, char *argv[])
{

    if (argc != 5)
    {
        cout<<"Usage: NOrMAL <config.txt> <forward_tags> <reverse_tags> <output_soft> <results>"<<endl;
        cout<<"       <config.txt>   - configuration file"<<endl;
        cout<<"       <forward_tags> - list of tags mapped on the forward strand"<<endl;
        cout<<"       <reverse_tags> - list of tags mapped on the reverse strand"<<endl;
        cout<<"       <results>      - output nucleosomes"<<endl;
        return 0;
    }
    Nucs RES;

    RES.Config( argv[1] ); //Read config file

    time_t before0, before,after ;
    before0 = time( NULL );
    std::cout<< "start reading input files" << std::endl;
    if( RES.Load_Data_X( argv[2]) && RES.Load_Data_Y( argv[3]) )
    {
        std::cout << "size of list X: " << RES.X.size() << std::endl;
        std::cout << "size of list Y: " << RES.Y.size() << std::endl;

        int start = std::min( RES.X.front(), RES.Y.front() );
        int finish = std::max( RES.X.back(), RES.Y.back() );
        int k = (finish - start)/ RES.NucSizePrior(); // #nucleosomes that are placed initially


        std::cout<<"Stage I"<<std::endl;
        RES.Init( start, finish, k); 
        
        bool flag = true;
        int count;
        int ind = 0;
        RES.PriorSet(0.001); //Specifying prior for nucleosome sizes, the lesser the value the stronger the prior
        while (flag) //Initial learning with strong prior for the nucleosome sizes to estimate number of nucleosomes
        {
            before = time( NULL);
            ind++;
            flag = false;
            count = 0;
            for (int i=0; i<7; i++)
                count += RES.EM_update();
            
            count += RES.Refine();
            if (count > 30 )
                flag = true;

            after = time( NULL );
            std::cout<< count << " clusters removed ";
            std::cout<< "cycle took " << after - before <<" sec"<< std::endl;
        }
        std::cout<<"Stage II relaxing prior"<<std::endl;
        flag = true;
        ind = 0;
        RES.PriorSet( RES.PriorEnd() ); //Relaxing prior nucleosome size to specified level to properly estimate nucleosome sizes
        while (flag)
        {
            before = time( NULL);
            ind++;
            flag = false;
            count = 0;
            for (int i=0; i<4; i++)
                count += RES.EM_update();
            
            count += RES.Refine();
            if (count > 30 )
                flag = true;

            after = time( NULL );
            std::cout<< count << " clusters removed ";
            std::cout<< "cycle took " << after - before <<" sec"<< std::endl;
        }
//        RES.Print( "dump.txt" ); //Uncomment this line if you want to dump nucleosome map before hardEM step
        std::cout<<"Stage III Hard EM"<<std::endl;
        std::cout<< "HardEM iteration removed" << RES.hardEM_update() << std::endl;
        after = time(NULL);
        std::cout<< RES.Param.size() << " clusters " << std::endl;
        RES.Print( argv[4]);
        std::cout<<" took overal time: "<< after - before0 <<" sec"<< std::endl;
    }
    return 0;
}
