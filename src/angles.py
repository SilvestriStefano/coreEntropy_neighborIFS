"""
Module that contains dynamics properties of a rational angle
"""
import numpy as np

from fractions import Fraction as Frac
from sympy import *


class Angle:
    def __init__(self,num,den):
        """       
        Parameters
        ----------
        num: int
            the numerator
        
        den: int
            the denominator

        """
        self.num = num
        self.den = den

    def orbit(self):
        """
        calculates the orbit of the angle num/den under the doubling map.

        Returns
        -------
        list
            the elements in the orbit of the angle
        
        """
        oList=[Frac(self.num,self.den)]
        while len(set(oList)) == len(oList):
            oList.append(oList[-1]*Frac(2,1)%1) #double it modulo 1
        return oList
    
    def period(self,orbit=None):
        """
        calculates the period of the angle num/den under the doubling map.

        Arguments
        ---------
        orbit: list
            the orbit of the angle

        Returns
        -------
        list
            the first entry is the index of where the period starts in the orbit 
            the second entry is the length of the period
        
        """
        orb = self.orbit() if orbit is None else orbit
        
        initPer = orb.index(orb[-1]) #where the period starts
        perLen = (len(orb)-initPer)-1 #length of the period in ks
        return [initPer,perLen]


    def ks_from_angle(self,orbit=None,period=None):
        """
        finds the kneading sequence associated to the angle num/den.

        Arguments
        ---------
        orbit: list
            the orbit of the angle
        period: list
            the index indicating where the period of the orbit starts and the length of it


        Returns
        -------
        str
            a string consisting of 1s and 0s
        
        """
        ks=''
        orb = self.orbit() if orbit is None else orbit
        perLen = self.period(orb)[1] if period is None else period[1]
        
        for num in orb:
            if num<Frac(self.num+self.den,2*self.den):
                ks+='1'
            else:
                ks+='0'

        if perLen==1:
            ks=ks[:-1]
            ks+='0'
        return ks
    
    def attr_itin_from_ks(self,ks,orbit=None,period=None):
        """
        finds the itinerary of 0 in the attractor given the kneading sequence 
        associated to the angle num/den.

        an entry of 1 in the ks corresponds to a change of sign in the itinerary,
        while an entry of 0 in the ks corresponds to the previous sign.

        Parameters
        ----------
        ks: str
            the kneading sequence

        Arguments
        ---------
        orbit: list
            the orbit of the angle
        period: list
            the index indicating where the period of the orbit starts and the length of it


        Returns
        -------
        str
            a string consisting of +s and -s
        
        """
        itin='+-' #we assume that ks starts with 1
        lastSign='-' #the current last symbol in the itinerary

        for i in range(1,len(ks)):
            vi=ks[i]
            if vi=='1':
                if lastSign=='+':
                    itin+='-'
                    lastSign='-'
                else:
                    itin+='+'
                    lastSign='+'
            else: #if vi=0
                if lastSign=='+':
                    itin+='+'
                    lastSign='+'
                else:
                    itin+='-'
                    lastSign='-'
        
        orb = self.orbit() if orbit is None else orbit
        initPer = self.period(orb)[0] if period is None else period[0]
        perLen = self.period(orb)[1] if period is None else period[1]

        if perLen==1: 
            #if the period of ks is 1, then the itin is complete, so remove the last char
            itin=itin[:-1]
        else:
            #otherwise check if the period of the itin is twice the one of the ks
            temp=''
            for i in range(initPer,len(ks)):
                vi=ks[i]
                if vi=='1':
                    if lastSign=='+':
                        temp+='-'
                        lastSign='-'
                    else:
                        temp+='+'
                        lastSign='+'
                else: #if vi=0
                    if lastSign=='+':
                        temp+='+'
                        lastSign='+'
                    else:
                        temp+='-'
                        lastSign='-'

            if temp!=itin[initPer+1:]:
                itin+=temp

        return itin

    def period_length_itin(self,itin,period):
        """
        finds the period length of the itinerary of 0 in the attractor.

        Parameters
        ----------
        itin: str
            the itinerary of 0

        Arguments
        ---------
        period: list
            the index indicating where the period of the orbit starts and the length of it

        Returns
        -------
        int
            the length of the period
        
        """

        initPer = period[0]
        perLen = period[1]
        perLenItin = perLen if perLen==1 else (len(itin)-initPer-1)

        return perLenItin

    def itin_to_rat(self, itin, period, pow_symb = '**'):
        """
        finds the numerator of the rational function associated to the itinerary of 0.
        
        Parameters
        ----------
        itin: str
            the itinerary of 0

        Arguments
        ---------
        period: list
            the index indicating where the period of the orbit starts and the length of it
        pow_symb: str
            the symbol for exponentiation. Default is **. Use ^ if you need to use 
            the result in Mathematica.

        Returns
        -------
        str
            a string containing a polynomial in x.
        
        """

        perLenItin = self.period_length_itin(itin,period)
        initPer = period[0]
        
        rat='('
        
        for pos in range(len(itin)):
            if ((perLenItin==1) & (pos==initPer)) :
                rat+=')*(1-x) +('+itin[pos]+'x'+pow_symb+str(pos)+')'
            elif pos==initPer+1:
                rat+=')*(1-x'+pow_symb+str(perLenItin)+') +('+itin[pos]+'x'+pow_symb+str(pos)
            elif pos==len(itin)-1:
                rat+=itin[pos]+'x'+pow_symb+str(pos)+')'
            else:
                rat+=itin[pos]+'x'+pow_symb+str(pos)
        
        return rat

    def assoc_lambda(self,rat_func):
        """
        find the roots of the associated rational function inside the disc of radius 2^(-1/2)

        Parameters
        ----------
        rat_func: str
            the numerator of the rational function
        
        Returns
        -------
        set
            the set of roots

        """

        rat = rat_func
        
        x = Symbol('x')
        f = Function('f')(x)

        f = simplify(rat)
        allroots = solveset(f)
        la = Intersection(allroots,ConditionSet(x,( np.abs(x)<=(1/np.sqrt(2)) ), S.Complexes ))
        try:
            lam = la.args[0]
            # raise Exception("possibly no roots")
        except:
            return la
        return lam