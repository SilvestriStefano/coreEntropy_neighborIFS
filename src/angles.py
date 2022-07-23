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
        self.frac = Frac(self.num,self.den)

    def __repr__(self):
        return f"Angle({self.num},{self.den})"
        
    def __str__(self):
        return f"The angle is {self.frac}"

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
        try: 
            orb = self.orbit_list
        except:
            # print("period()->orbit was not called")
            self.orbit()
            orb = self.orbit_list
            # print("----------period()------------")
        
        self.start_index_per = orb.index(orb[-1]) #where the period starts in the orbit
        self.per_len = (len(orb)-self.start_index_per)-1 #length of the period in ks


    def ks_from_angle(self):
        """
        finds the kneading sequence associated to the angle num/den.
        
        Use the attribute `ks` to get the result
        
        """
        self.ks=''
        try: 
            orb = self.orbit_list
        except:
            # print("  ks_from_angle()->orbit was not called")
            self.orbit()
            orb = self.orbit_list
            # print("  ---------ks_from_angle()1------------")

        for numb in orb:
            if Frac(self.num,2*self.den) < numb < Frac(self.num+self.den,2*self.den):
                self.ks+='1'
            elif numb == Frac(self.num,2*self.den) or numb == Frac(self.num+self.den,2*self.den):
                self.ks+='*'
            else:
                self.ks+='0'
        
        try: 
            period_length = self.per_len
        except:
            # print("  ks_from_angle()->period was not called")
            self.period()
            period_length = self.per_len
            # print("  ----------ks_from_angle()2------------")

        self.ks=self.ks[:-1]
        # if period_length==1:
        #     self.ks=self.ks[:-1]
        #     self.ks+='0'
        # else:
        #     self.ks=self.ks[:-1]

    
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

        try:
            kn_seq = self.ks
        except:
            # print("    attr_itin_from_ks()->ks_from_angle was not called")
            self.ks_from_angle()
            kn_seq = self.ks
            # print("    ---------------attr_itin_from_ks()1--------------")


        for vi in kn_seq[1:]:
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

        try: 
            period_length = self.per_len
            starting_index = self.start_index_per
        except:
            # print("    attr_itin_from_ks()->either orbit and-or period was not called")
            self.orbit()
            self.period()
            period_length = self.per_len    
            starting_index = self.start_index_per
            # print("    -------------------attr_itin_from_ks()2-----------------------")

        if period_length==1:
            #if the period of ks is 1, then the itin is complete, so remove the last char
            self.itin=self.itin[:-1]
        else:
            #otherwise check if the period of the itin is twice the one of the ks
            temp=''
            for vi in kn_seq[starting_index:]:
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

            if temp!=self.itin[starting_index+1:]:
                self.itin+=temp


    def period_length_itin(self):
        """
        finds the period length of the itinerary of 0 in the attractor.
        Use the attribute `perLenItin` to get the value.
        """
        try:
            it = self.itin
        except:
            # print("      period_length_itin()->attr_itin_from_ks was not called")
            self.attr_itin_from_ks()
            it = self.itin
            # print("      -----------period_length_itin()-----------------------")

        self.itin_per_len = self.per_len if self.per_len==1 else (len(it)-self.start_index_per-1)


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
        try:
            it = self.itin
        except:
            # print("        itin_to_rat()->attr_itin_from_ks was not called")
            self.attr_itin_from_ks()
            it = self.itin
            # print("        ---------------itin_to_rat()1------------------")
        try:
            it_period = self.itin_per_len
        except:
            # print("        itin_to_rat()->period_length_itin was not called")
            self.period_length_itin()
            it_period = self.itin_per_len
            # print("        ----------------itin_to_rat()2------------------")

        self.rat_func = '('
        
        for index, sign in enumerate(it):
            if ((it_period==1) & (index==self.start_index_per)) :
                self.rat_func += ')*(1-x) +('+sign+'x'+pow_symb+str(index)+')'
            elif (index==self.start_index_per+1):
                self.rat_func += ')*(1-x'+pow_symb+str(it_period)+') +('+sign+'x'+pow_symb+str(index)
            elif index==len(it)-1:
                self.rat_func += sign+'x'+pow_symb+str(index)+')'
            else:
                self.rat_func += sign+'x'+pow_symb+str(index)
        

    def assoc_lambda(self):
        """
        find the roots of the associated rational function inside the disc of radius 2^(-1/2)
        
        Use the attribute `lam` to get the result
        """
        try:
            num_poly = self.rat_func.replace('^','**')
        except:
            # print("          assoc_lambda()->itin_to_rat was not called")
            self.itin_to_rat()
            num_poly = self.rat_func
            # print("          --------------assoc_lambda()--------------")

        x = Symbol('x')
        f = Function('f')(x)

        f = simplify(num_poly)
        allroots = solveset(f)
        la = Intersection(allroots,ConditionSet(x,( np.abs(x)<=(1/np.sqrt(2)) ), S.Complexes ))
        try:
            self.lam = la.args[0]
            # raise Exception("possibly no roots")
        except:
            self.lam = la
