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

        Use the attribute `orbit_list` to get the result
        
        """
        self.orbit_list=[Frac(self.num,self.den)]
        while len(set(self.orbit_list)) == len(self.orbit_list):
            self.orbit_list.append(self.orbit_list[-1]*Frac(2,1)%1) #double it modulo 1
    
    def period(self):
        """
        calculates the period of the angle num/den under the doubling map and 
        the index in the orbit list of where it starts.

        Use the attribute `start_index_per` to get the index.
        Use the attribute `perLen` to get the length of the period.
        
        """
        orb = self.orbit_list
        
        self.start_index_per = orb.index(orb[-1]) #where the period starts in the orbit
        self.per_len = (len(orb)-self.start_index_per)-1 #length of the period in ks


    def ks_from_angle(self):
        """
        finds the kneading sequence associated to the angle num/den.
        
        Use the attribute `ks` to get the result
        
        """
        self.ks=''
        orb = self.orbit_list
        
        for numb in orb:
            if numb<Frac(self.num+self.den,2*self.den):
                self.ks+='1'
            else:
                self.ks+='0'

        if self.per_len==1:
            self.ks=self.ks[:-1]
            self.ks+='0'

    
    def attr_itin_from_ks(self):
        """
        finds the itinerary of 0 in the attractor given the kneading sequence 
        associated to the angle num/den.

        an entry of 1 in the ks corresponds to a change of sign in the itinerary,
        while an entry of 0 in the ks corresponds to the previous sign.

        Use the attribute `itin` to get the result
        
        """
        self.itin='+-' #we assume that ks starts with 1
        lastSign='-' #the current last symbol in the itinerary

        for i in range(1,len(self.ks)):
            vi=self.ks[i]
            if vi=='1':
                if lastSign=='+':
                    self.itin+='-'
                    lastSign='-'
                else:
                    self.itin+='+'
                    lastSign='+'
            else: #if vi=0
                if lastSign=='+':
                    self.itin+='+'
                    lastSign='+'
                else:
                    self.itin+='-'
                    lastSign='-'
        
        orb = self.orbit_list
        
        if self.per_len==1: 
            #if the period of ks is 1, then the itin is complete, so remove the last char
            self.itin=self.itin[:-1]
        else:
            #otherwise check if the period of the itin is twice the one of the ks
            temp=''
            for i in range(self.start_index_per,len(self.ks)):
                vi=self.ks[i]
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

            if temp!=self.itin[self.start_index_per+1:]:
                self.itin+=temp


    def period_length_itin(self):
        """
        finds the period length of the itinerary of 0 in the attractor.
        Use the attribute `perLenItin` to get the value.
        """
        self.itin_per_len = self.per_len if self.per_len==1 else (len(self.itin)-self.start_index_per-1)


    def itin_to_rat(self, pow_symb = '**'):
        """
        finds the numerator of the rational function associated to the itinerary of 0.

        Use the attribute `rat_func` to get the result

        Arguments
        ---------
        pow_symb: str
            the symbol for exponentiation. Default is **. Use ^ if you need to use 
            the result in Mathematica.

        """

        self.rat_func = '('
        
        for pos in range(len(self.itin)):
            if ((self.itin_per_len==1) & (pos==self.start_index_per)) :
                self.rat_func += ')*(1-x) +('+self.itin[pos]+'x'+pow_symb+str(pos)+')'
            elif pos==self.start_index_per+1:
                self.rat_func += ')*(1-x'+pow_symb+str(self.itin_per_len)+') +('+self.itin[pos]+'x'+pow_symb+str(pos)
            elif pos==len(self.itin)-1:
                self.rat_func += self.itin[pos]+'x'+pow_symb+str(pos)+')'
            else:
                self.rat_func += self.itin[pos]+'x'+pow_symb+str(pos)
        

    def assoc_lambda(self):
        """
        find the roots of the associated rational function inside the disc of radius 2^(-1/2)
        
        Use the attribute `lam` to get the result
        """

        x = Symbol('x')
        f = Function('f')(x)

        f = simplify(self.rat_func)
        allroots = solveset(f)
        la = Intersection(allroots,ConditionSet(x,( np.abs(x)<=(1/np.sqrt(2)) ), S.Complexes ))
        try:
            self.lam = la.args[0]
            # raise Exception("possibly no roots")
        except:
            self.lam = la