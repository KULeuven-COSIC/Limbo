# Limbo
[Limbo](https://eprint.iacr.org/2021/215) implementation


## Setup
The config.h file lets you define some constants that slightly change the behaviour of the program.

Security levels and other parameters are defined in limbo_instance.cpp
You have to specify yours according to the circuit you want to prove.
In particular, if you do not give the correct number of mult gates, 
as well as input and output wires this will result in an error.

By default, the parameter set 0 works for the SHA256 circuit.

To run limbo:
```bash
mkdir build
cd build
cmake ..
make 
wget https://homes.esat.kuleuven.be/~nsmart/MPC/sha256.txt
# benchmarks
./bench sha256.txt 0
```


More circuits in the correct format are available [here](https://homes.esat.kuleuven.be/~nsmart/MPC/)


## Parameters
As for all MPCitH based system, there is a trade-off between computational complexity and proof size.
The more you increase the number of in the head parties, the lower the number of repetitions needs to be.
Thus reducing the proof size while increasing the computational complexity.


In the provided example for SHA256 (parameter set 0 in limbo_instance.cpp) we use 64 parties and 29 repetitions to achieve 128 bit of security.
If the goal was to otpimize for speed, we could for example choose 8 parties and 48 repetitions.


In [our paper](https://eprint.iacr.org/2021/215) we provide parameters set for different circuit size, different trade-off
and different target security.

## Acknowledgements

Some parts of the code were based on the [optimized Picnic implementation](https://github.com/IAIK/Picnic) and on the [banquet implementation](https://github.com/dkales/banquet).
