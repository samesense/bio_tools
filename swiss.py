import utils_swiss

seq = 'MRIILLGAPGAGKGTQAQFIMAKFGIPQISTGDMLRAAIKAGTELGKQAK'
[percent, pdb] = utils_swiss.getSwissModelPage(seq)
print percent
