#include "Scoring.hpp"
#include <iostream>
#include "RealOptions.hpp"

int main(int argc, char *argv[])
{
        try
        {
        	RealOptions opts(argc,argv);
                Scoring const scoring(opts.similarity, opts.gc, opts.trans, opts.err, opts.gcmut_bias);
                std::cout << scoring;
                return EXIT_SUCCESS;
        }
        catch(std::exception const & ex)
        {
                std::cerr << ex.what() << std::endl;
                return EXIT_FAILURE;
        }
}
