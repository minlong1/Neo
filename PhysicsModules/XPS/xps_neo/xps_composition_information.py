#File of all literature information into the BE, sigma/gamma, and spin-orbitting splitting ranges for specific elemental compositions
#Could also just move all this information into the periodic table file
class ElementalComposition:
    def __init__(self,element,photoelectronLine, BE_guess):
        elementalComps = ['C', 'SiC', 'CO2'] #This is for testing
        self.elementalComps = elementalComps #--> Make input array for Fitting parameters tab that appends this info
        #Data Presented as C, SiC, CO2 etc.
        self.element = element
        self.photoelectronLine = photoelectronLine
        #self.BE_lit, self.is_singlet, self.so_split, self.branching_ratio = self.getParams(self.element,photoelectronLine)
        self.BE_guess = BE_guess

    def getParamRanges(self):
        if self.element == 'C': #Should I make seperate function calls for each element or just have a lot of if statements?

            for i in self.elementalComps:
                if i == 'C':
                    BE_guess[i] = 284.5 #We want to append these values to an array for each peak. Call BE_guess from xps.py?
                elif i == 'SiC':
                    BE_guess[i] = 282.0
