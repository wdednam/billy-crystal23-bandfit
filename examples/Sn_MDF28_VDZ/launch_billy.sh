#PBS -o output.txt
#PBS -e error.txt
#PBS -q workq
#PBS -l software=test
#PBS -l select=1
#PBS -l walltime=120:00:00
#PBS -l place=pack
cd $PBS_O_WORKDIR
source /home/dednaw/CRYSTAL23/cry23.bashrc
# Launch program
billy -f -n 3 -np 1 -cost ref.band.npy -wE 0.2 -lambda 0.005 -cap 0.75 -erange -10 10 Sn 15
