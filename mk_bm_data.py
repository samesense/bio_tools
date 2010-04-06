import random

with open('junk', 'w') as f:
    for x in xrange(5000):
        f.write('gene' + str(x) + '\t' + str(random.gauss(.1, .05)) + '\n')
        f.write('gene' + str(x) + '\t' + str(random.gauss(.3, .05)) + '\n')
#os.system('python ../Figure/bm_trainer.py junk > junk.out')
