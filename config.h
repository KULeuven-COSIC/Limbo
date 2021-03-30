//Can have up to 64*PARTY64 parties
//(tested only for PARTY64=1 or 2)
#define PARTY64 1

//Defines the set of parameters to use
//(for now we only have 40 or 128)
#define SEC 128


//Defines if we want to use same challenge
//(comment out for indep challenges)
//#define SAME_CHALLENGE

//Defines whether I want timers inside the code
//(comment out to remove)
//#define IN_CODE_TIMERS

//Defines whether I want very detailed timers
//(comment out to remove)
//#define DETAILED_IN_CODE_TIMERS

//Defines whether we perform some sanity checks
//(in particular, prover checks he ends up with inner product)
//(comment out to not do those)
#define SANITY_CHECK
