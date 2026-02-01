#PBS -o output.txt
#PBS -e error.txt
#PBS -q workq
#PBS -l software=test
#PBS -l select=1
#PBS -l walltime=120:00:00
#PBS -l place=pack
cd $PBS_O_WORKDIR
source $HOME/CRYSTAL23/cry23.bashrc
# Launch program
billy -f -n 3 -np 1 -cost ref.band.npy -wE 0.15 -lambda 0.01 -cap 0.55 -erange -10.0 10.0 Al 15
