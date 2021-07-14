class TimeHandler( object ):
    def __init__(self, timeStep, finalTime, unit="Second", initialTime=0.0 ):
        self.unit = unit
        self.__timeStep = self.getTime(timeStep)
        self.__finalTime = self.getTime(finalTime)
        self.__initialTime = self.getTime(initialTime)
        self.__currentTime = self.getTime(initialTime)
        if self.__timeStep > self.__finalTime:
            raise NameError("Time step size cannot be larger than final time.")

    def advanceTime(self):
        self.__currentTime += self.__timeStep

    def getTimeStep(self):
        return self.__timeStep

    def setTimeStep(self, timeStep):
        self.__timeStep = timeStep

    def getCurrentTime(self):
        return self.__currentTime

    def getFinalTime(self):
        return self.__finalTime

    def getInitialTime(self):
        return self.__initialTime

    def isFinalTimeReached( self ):
        if self.__currentTime >= self.__finalTime:
            return False
        else:
            return True

    def printCurrentTime(self):
        print(self.__currentTime)

    def getTime(self, time):
        if self.unit == "Second":       conversion = 1
        elif self.unit == "Minute":     conversion = 60
        elif self.unit == "Hour":       conversion = 60*60
        elif self.unit == "Day":        conversion = 24*60*60
        elif self.unit == "Week":       conversion = 7*24*60*60
        elif self.unit == "Year":       conversion = 365*24*60*60
        else:                           raise NameError("Time unit %s no supported."%self.unit)
        return time*conversion

class IterativeCycleController(object):
    def __init__(self, maxIte, maxTol):
        self.maxIte = maxIte
        self.maxTol = maxTol
        self.__keepCycle = True
        self.iteNumber = 0

    def keepCycling(self):
        return self.__keepCycle

    def printKeepCycle(self):
        print(self.__keepCycle)

    def reset(self):
        self.__keepCycle = True
        self.iteNumber = 0

    def execute(self, error):
        if error < self.maxTol or self.iteNumber > self.maxIte:
            if self.iteNumber > 1:
                self.__keepCycle = False
        self.iteNumber += 1
