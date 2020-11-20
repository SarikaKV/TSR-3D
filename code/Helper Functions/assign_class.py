

def thetaClass_(binBoundaries, Theta, type):
	for i in binBoundaries:
		if Theta < binBoundaries[0]:
			if type == 0:
				print 'out of index'
			else:
				print binBoundaries.index(binBoundaries[0])+1
			break
		if (Theta < i) :
			if type ==0:
				print 'At '+ str(binBoundaries.index(i) )
			else:
				print 'At '+ str(binBoundaries.index(i) + 1)
			break
	if Theta >= binBoundaries[-1]:
		if type ==0:
			if Theta == binBoundaries[-1]:
				print 'At '+ str(binBoundaries.index(binBoundaries[-1]))
		else :
			print 'At '+ str(binBoundaries.index(binBoundaries[-1]) +2)
	
	


thetaBounds = [ 0.1 , 26.68, 39.82, 51.1, 61.4, 71.16, 80.64, 90]
dTheta = len(thetaBounds)
print dTheta -1
distBounds =  [2.81, 7.00, 9.00, 11.00, 14.00, 17.4, 24.16, 30.19, 36.37, 44.78, 175.52]
dLen = len(distBounds) + 1
print dLen
Theta = 90
thetaClass_(thetaBounds,Theta,0)
thetaClass_(distBounds,175.52,1)