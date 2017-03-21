from joblib import Parallel, delayed
from itertools import product
from os.path import isfile
from os import system
from sys import argv
def runHM(seq_rate_qp_fskip,cfg,numFrames):
	[seq_rate,qp,fskip] = seq_rate_qp_fskip
	[seq,rate] = seq_rate
	fskip_idx = fskip*rate
	

	debugPath = "./debug-out/%s_qp%s_%dfr_fskip%d.csv" % (seq, qp, numFrames,fskip_idx)
	if not isfile(debugPath):
		cmd = './hm-debug -c ../cfg/%s.cfg -c ~/hm-cfgs/server/%s.cfg -q %s -f %d --FrameSkip=%d --DebugFile=%s' % (cfg,seq, qp, numFrames,fskip_idx, debugPath)
		print cmd
		system(cmd)
	


sequences_A = [['Traffic',30], ['SteamLocomotiveTrain_10bit',60], ['NebutaFestival_10bit',60], ['PeopleOnStreet',30]]
sequences_B = [['BQTerrace',60], ['ParkScene',24], ['BasketballDrive',50], ['Cactus',50], ['Kimono',24],['Tennis',24 ]]
sequences_720p = [['Johnny',60], ['FourPeople',60], ['SlideShow',20],['ChinaSpeed',30], ['SlideEditing',30]]
sequences_C = [['BasketballDrill',50], ['RaceHorsesC',32], ['PartyScene',50], ['BQMall',60]]
sequences_D = [['BasketballPass',50],  ['RaceHorses',32],  ['BlowingBubbles',50], ['BQSquare',60], ['FlowerVase', 32]]
#sequences = [['RaceHorses',30]]
sequences = sequences_A + sequences_B
#sequences = sequences_C + sequences_D
try:
	cfg = argv[1]
	numFrames = int(argv[2])
	n_cores = int(argv[3])
except:
	print 'python %s CFG NUM_FRAMES N_CORES' % argv[0]
	exit()
qps = ['22','27','32','37','42']
#qps = ['27','32']
secs = [0,1,2,3,4]
secs = [0,1,2,3,4]

system('mkdir -p ./debug-out')
#sequences = [['BQSquare',60]]
#secs = [0]

Parallel(n_jobs=n_cores)(delayed(runHM)(seq_rate_qp_fskip,cfg,numFrames) for seq_rate_qp_fskip in product(sequences, qps,secs))
