"""
Module that contains dynamics properties of a rational angle
"""
from os import path

import numpy as np
from fractions import Fraction as Frac
from sympy import Symbol, S, Function, Intersection, ConditionSet, Abs, simplify, solveset
from sympy.core import Add

import logging
import logging.config

src_dir, _ = path.split(path.abspath(__file__))
log_conf_path = path.join(path.dirname(src_dir),'log/logging.conf')
logging.config.fileConfig(log_conf_path)

# create logger
logger = logging.getLogger("default")


class Angle:
    def __init__(self,num:int=None,den:int=None,*,th:str=None):
        """
        Choose between initializing with two integers or a string

        Parameters
        ----------
        num: int
            the numerator
        den: int
            the denominator
        th: str
            The rational angle in string format
        """

        if num is not None and den is not None:
            self.num = num
            self.den = den
        if th is not None:
            theta = th.split("/")
            try:
                self.num = int(theta[0])
                self.den = int(theta[1])
            except ValueError:
                raise ValueError("Please provide a valid angle string (e.g. '3/4')") from None
            except IndexError:
                raise ValueError("Please provide a valid angle string (e.g. '3/4')") from None

        self.frac = Frac(self.num,self.den)

    def __repr__(self):
        return f"Angle({self.num},{self.den})"
        
    def __str__(self):
        return f"{self.num}/{self.den}"

    def orbit(self)->list:
        """
        calculates the orbit of the angle num/den under the doubling map.

        Return
        ------
        orbit_list: list
            It is a list of Frac elements.
        
        """
        if getattr(self,"orbit_list",None): # avoid recalculating if it already exists
            return self.orbit_list
        
        self.orbit_list=[Frac(self.num,self.den)]
        while len(set(self.orbit_list)) == len(self.orbit_list):
            self.orbit_list.append(self.orbit_list[-1]*Frac(2,1)%1) #double it modulo 1
        return self.orbit_list
    
    def period(self)->tuple:
        """
        calculates the period of the angle num/den under the doubling map and 
        the index in the orbit list of where it starts.

        Returns
        -------
        tuple of int
            The first entry of the tuple is `per_len`, the length of the period.
            The second entry of the tuple is `start_index_per`, the index of where the periodic part of the orbit starts.
        
        """

        if getattr(self,"per_len",None) and getattr(self,"start_index_per",None): # avoid recalculating if it already exists
            return (self.per_len,self.start_index_per)
        
        orb = self.orbit()
        
        self.start_index_per = orb.index(orb[-1]) #where the period starts in the orbit
        self.per_len = (len(orb)-self.start_index_per)-1 #length of the period in ks
        return (self.per_len,self.start_index_per)

    def ks_from_angle(self)->str:
        """
        finds the kneading sequence associated to the angle num/den.
        
        Return
        ------
        ks: string
            The kneading sequence.
        
        """

        if getattr(self,"ks",None): # avoid recalculating if it already exists
            return self.ks
        
        self.ks=''
        orb = self.orbit()

        for numb in orb:
            if Frac(self.num,2*self.den) < numb < Frac(self.num+self.den,2*self.den):
                self.ks+='1'
            elif numb == Frac(self.num,2*self.den) or numb == Frac(self.num+self.den,2*self.den):
                self.ks+='*'
            else:
                self.ks+='0'

        self.ks=self.ks[:-1]
        return self.ks
    
    def attr_itin_from_ks(self)->str:
        """
        finds the itinerary of 0 in the attractor given the kneading sequence 
        associated to the angle num/den.

        an entry of 1 in the ks corresponds to a change of sign in the itinerary,
        while an entry of 0 in the ks corresponds to the previous sign.

        Return
        ------
        itin: string
            The itinerary of 0 in the attractor.
        
        """

        if getattr(self,"itin",None): # avoid recalculating if it already exists
            return self.itin
        
        self.itin='+-' #we assume that ks starts with 1
        lastSign='-' #the current last symbol in the itinerary

        kn_seq = self.ks_from_angle()

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

        period_length, starting_index = self.period()

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
        return self.itin

    def period_length_itin(self)->int:
        """
        finds the period length of the itinerary of 0 in the attractor.
        
        Return
        ------
        itin_per_len: int
            The length of the periodic part of the itinerary.
        """

        if getattr(self,"itin_per_len",None): # avoid recalculating if it already exists
            return self.itin_per_len
        
        it = self.attr_itin_from_ks()

        self.itin_per_len = self.per_len if self.per_len==1 else (len(it)-self.start_index_per-1)
        return self.itin_per_len

    def itin_to_rat(self, *, pow_symb:str = '**')->str:
        """
        finds the numerator of the rational function associated to the itinerary of 0.


        Parameters
        ----------
        pow_symb: str
            the symbol for exponentiation. Default is '**'. Use '^' if you need to use 
            the result in Mathematica.

        Return
        ------
        rat_func: string
            The numerator (polynomial) of the rational function which vanishes at lambda.

        """

        if getattr(self,"rat_func",None): # avoid recalculating if it already exists
            symb_split = self.rat_func.split(pow_symb)
            if len(symb_split)>1: 
                return self.rat_func
            else: 
                old_symb = "**" if pow_symb =="^" else "^"
                self.rat_func = self.rat_func.replace(old_symb,pow_symb)
                return self.rat_func
        
        it = self.attr_itin_from_ks()
        it_period = self.period_length_itin()

        self.rat_func = '('
        
        for index, sign in enumerate(it):
            if ((it_period==1) & (index==self.start_index_per)) :
                if self.rat_func=='(': # special case when angle is 0/1
                    self.rat_func = '1'
                else:
                    self.rat_func += ')*(1-x) +('+sign+'x'+pow_symb+str(index)+')'
            elif (index==self.start_index_per+1):
                self.rat_func += ')*(1-x'+pow_symb+str(it_period)+') +('+sign+'x'+pow_symb+str(index)
            elif index==len(it)-1:
                self.rat_func += sign+'x'+pow_symb+str(index)+')'
            else:
                self.rat_func += sign+'x'+pow_symb+str(index)
        return self.rat_func

    def assoc_lambda(self)->Add:
        """
        find the roots of the associated rational function inside the disc of radius 2^(-0.5)+10^(-14).
        If none is found it returns 0.+0.j.
        
        Return
        ------
        lam: sympy.core.add.Add
            The complex number *associated* to the Misiurewicz parameter.
            To get the value simply cast its type to complex: `complex(lam)`.
        """

        if getattr(self,"lam",None): # avoid recalculating if it already exists
            return self.lam
        
        num_poly = self.itin_to_rat()

        x = Symbol('x')
        f = Function('f')(x)

        f = simplify(num_poly)
        allroots = solveset(f)
        la = Intersection(allroots, ConditionSet(x,( Abs(x)<=(1/np.sqrt(2)+1e-14) ), S.Complexes))
        try:
            self.lam = la.args[0]
            # raise Exception("possibly no roots")
        except IndexError as e:
            logger.error(f"{self}; Could not find any viable solution inside the disk of radius 2^(-0.5)")
            self.lam = 0.+0.j
        return self.lam
