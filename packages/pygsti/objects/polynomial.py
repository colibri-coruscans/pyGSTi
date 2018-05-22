""" Defines the Polynomial class """
from __future__ import division, print_function, absolute_import, unicode_literals
#*****************************************************************
#    pyGSTi 0.9:  Copyright 2015 Sandia Corporation
#    This Software is released under the GPL license detailed
#    in the file "license.txt" in the top-level pyGSTi directory
#*****************************************************************

import numpy as _np
try:
    from . import fastreplib as replib
except ImportError:
    from . import replib



class Polynomial(dict):
    """ Encapsulates a polynomial """

    @classmethod
    def fromrep(cls, rep):
        """ TODO docstring """
        max_num_vars = rep.max_num_vars  # one of the few/only cases where a rep
        max_order = rep.max_order        # needs to expose some python properties

        def int_to_vinds(indx):
            ret = []
            while indx != 0:
                nxt = indx // (max_num_vars+1)
                i = indx - nxt*(max_num_vars+1)
                ret.append(i-1)
                indx = nxt
            assert(len(ret) <= max_order)
            return tuple(sorted(ret))
        
        tup_coeff_dict = { int_to_vinds(k): val for k,val in rep.coeffs.items() }
        return cls(tup_coeff_dict)

    
    def __init__(self, coeffs=None):
        """ TODO: docstring - coeffs is a dict of coefficients w/keys == tuples
             of integer variable indices.  E.g. (1,1) means "variable1 squared"
        """
        super(Polynomial,self).__init__()
        if coeffs is not None:
            self.update(coeffs)
            
    def deriv(self, wrtParam):
        dcoeffs = {}
        for ivar, coeff in self.items():
            cnt = float(ivar.count(wrtParam))
            if cnt > 0:
                l = list(ivar)
                del l[l.index(wrtParam)]
                dcoeffs[ tuple(l) ] = cnt * coeff

        return Polynomial(dcoeffs)

    def get_max_order(self):
        return max([len(k) for k in self.keys()])

    def evaluate(self, variable_values):
        """ TODO: docstring -- and make this function smarter (Russian peasant) """
        ret = 0
        for ivar,coeff in self.items():
            ret += coeff * _np.product( [variable_values[i] for i in ivar] )
        return ret

    def compact(self):
        """ TODO docstring Returns compact representation of (vtape, ctape) 1D nupy arrays """
        iscomplex = any([ abs(_np.imag(x)) > 1e-12 for x in self.values() ])
        nTerms = len(self)
        nVarIndices = sum(map(len,self.keys()))
        vtape = _np.empty(1 + nTerms + nVarIndices, 'i') # "variable" tape
        ctape = _np.empty(nTerms, complex if iscomplex else 'd') # "coefficient tape"

        i = 0
        vtape[i] = nTerms; i+=1
        for iTerm,k in enumerate(sorted(self.keys())):
            l = len(k)
            ctape[iTerm] = self[k] if iscomplex else _np.real(self[k])
            vtape[i] = l; i += 1
            vtape[i:i+l] = k; i += l
        assert(i == len(vtape)), "Logic Error!"
        return vtape, ctape

    def copy(self):
        return Polynomial(self)

    def map_indices(self, mapfn):
        """ TODO: docstring - mapfn should map old->new variable-index-tuples """
        new_items = { mapfn(k): v for k,v in self.items() }
        self.clear()
        self.update(new_items)

    def mult(self,x):
        """ Does self * x where x is a polynomial """
        newpoly = Polynomial()
        for k1,v1 in self.items():
            for k2,v2 in x.items():
                k = tuple(sorted(k1+k2))
                if k in newpoly: newpoly[k] += v1*v2
                else: newpoly[k] = v1*v2
        return newpoly

    def scale(self, x):
        # assume a scalar that can multiply values
        for k in tuple(self.keys()): # I think the tuple() might speed things up (why?)
            self[k] *= x

    def scalar_mult(self, x):
        newpoly = self.copy()
        newpoly.scale(x)
        return newpoly

    def __str__(self):
        def fmt(x):
            if abs(_np.imag(x)) > 1e-6:
                if abs(_np.real(x)) > 1e-6: return "(%.3f+%.3fj)" % (x.real, x.imag)
                else: return "(%.3fj)" % x.imag
            else: return "%.3f" % x.real
            
        termstrs = []
        sorted_keys = sorted(list(self.keys()))
        for k in sorted_keys:
            varstr = ""; last_i = None; n=0
            for i in sorted(k):
                if i == last_i: n += 1
                elif last_i is not None:
                    varstr += "x%d%s" % (last_i, ("^%d" % n) if n > 1 else "")
                last_i = i
            if last_i is not None:
                varstr += "x%d%s" % (last_i, ("^%d" % n) if n > 1 else "")
            #print("DB: k = ",k, " varstr = ",varstr)
            if abs(self[k]) > 1e-4:
                termstrs.append( "%s%s" % (fmt(self[k]), varstr) )
        if len(termstrs) > 0:
            return " + ".join(termstrs)
        else: return "0"

    def __repr__(self):
        return "Poly[ " + str(self) + " ]"

    def __add__(self,x):
        newpoly = self.copy()
        if isinstance(x, Polynomial):
            for k,v in x.items():
                if k in newpoly: newpoly[k] += v
                else: newpoly[k] = v
        else: # assume a scalar that can be added to values
            for k in newpoly:
                newpoly[k] += x
        return newpoly

    def __iadd__(self,x):
        """ Does self += x more efficiently """
        if isinstance(x, Polynomial):
            for k,v in x.items():
                try:
                    self[k] += v
                except KeyError:
                    self[k] = v
        else: # assume a scalar that can be added to values
            for k in self:
                self[k] += x
        return self

    def __mul__(self,x):
        #if isinstance(x, Polynomial):
        #    newpoly = Polynomial()
        #    for k1,v1 in self.items():
        #        for k2,v2 in x.items():
        #            k = tuple(sorted(k1+k2))
        #            if k in newpoly: newpoly[k] += v1*v2
        #            else: newpoly[k] = v1*v2
        #else:
        #    # assume a scalar that can multiply values
        #    newpoly = self.copy()
        #    for k in newpoly:
        #        newpoly[k] *= x
        #return newpoly
        if isinstance(x, Polynomial):
            return self.mult(x)
        else: # assume a scalar that can multiply values
            return self.scale(x)

    def __rmul__(self, x):
        return self.__mul__(x)

    def __pow__(self,n):
        ret = Polynomial({(): 1.0}) # max_order updated by mults below
        cur = self
        for i in range(int(np.floor(np.log2(n)))+1):
            rem = n % 2 #gets least significant bit (i-th) of n
            if rem == 1: ret *= cur # add current power of x (2^i) if needed  
            cur = cur*cur # current power *= 2
            n //= 2 # shift bits of n right 
        return ret

    def __copy__(self):
        return self.copy()

    def torep(self, max_order=None, max_num_vars=None):
        """ TODO: docstring """
        # Set max_order (determines based on coeffs if necessary)
        default_max_order = 0 if len(self) == 0 else \
                            max([ len(k) for k in self.keys()])
        if max_order is None:
            max_order = default_max_order
        else:
            assert(default_max_order <= max_order)

        # Set max_num_vars (determines based on coeffs if necessary)
        default_max_vars = 0 if len(self) == 0 else \
                           max([ (max(k)+1 if k else 0) for k in self.keys()])
        if max_num_vars is None:
            max_num_vars = default_max_vars
        else:
            assert(default_max_vars <= max_num_vars)

        #new.max_order = max_order            
        #new.max_num_vars = max_num_vars
        def vinds_to_int(vinds):
            """ Convert tuple index of ints to single int given max_order,max_numvars """
            assert(len(vinds) <= max_order), "max_order is too low!"
            ret = 0; m = 1
            for i in vinds: # last tuple index is most significant                                                                                                          
                assert(i < max_num_vars), "Variable index exceed maximum!"
                ret += (i+1)*m
                m *= max_num_vars+1
            return ret
        
        int_coeffs = { vinds_to_int(k): v for k,v in self.items() }
        return replib.PolyRep(int_coeffs, max_order, max_num_vars)
    

#OLD: TODO REMOVE
#class SLOWPolynomial(object):
#    """ Encapsulates a polynomial """
#    def __init__(self, coeffs=None,is_complex=True):
#        """ TODO: docstring - coeffs is a dict of coefficients w/keys == tuples
#             of integer variable indices.  E.g. (1,1) means "variable1 squared"
#        """
#        if coeffs is None:
#            self.coeffs = _np.zeros(0,complex if is_complex else 'd')
#            self.inds = []
#        else:
#            self.inds = sorted(list(coeffs.keys()))
#            self.coeffs = _np.array([coeffs[k] for k in self.inds],complex if is_complex else 'd')
#            
#    def deriv(self, wrtParam):
#        dcoeffs = {}
#        for ivar,coeff in zip(self.inds,self.coeffs):
#            cnt = float(ivar.count(wrtParam))
#            if cnt > 0:
#                l = list(ivar)
#                del l[ivar.index(wrtParam)]
#                dcoeffs[ tuple(l) ] = cnt * coeff
#
#        return Polynomial(dcoeffs) # returns another polynomial
#
#    def evaluate(self, variable_values):
#        """ TODO: docstring -- and make this function smarter (Russian peasant) """
#        ret = 0
#        for ivar,coeff in zip(self.inds,self.coeffs):
#            ret += coeff * _np.product( [variable_values[i] for i in ivar] )
#        return ret
#
#    def compact(self):
#        """ TODO docstring Returns compact representation of (vtape, ctape) 1D nupy arrays """
#        iscomplex = bool(_np.linalg.norm(_np.imag(self.coeffs)) > 1e-12)
#        nTerms = len(self.inds)
#        nVarIndices = sum(map(len,self.inds))
#        vtape = _np.empty(1 + nTerms + nVarIndices, 'i') # "variable" tape
#        ctape = _np.empty(nTerms, complex if iscomplex else 'd') # "coefficient tape"
#
#        i = 0
#        vtape[i] = nTerms; i+=1
#        for iTerm,k in enumerate(self.inds):
#            l = len(k)
#            ctape[iTerm] = self.coeffs[iTerm] if iscomplex else _np.real(self.coeffs[iTerm])
#            vtape[i] = l; i += 1
#            vtape[i:i+l] = k; i += l
#        assert(i == len(vtape)), "Logic Error!"
#        return vtape, ctape
#
#    def copy(self):
#        cpy = Polynomial()
#        cpy.coeffs = self.coeffs.copy()
#        cpy.inds = self.inds[:]
#        return cpy
#
#    def map_indices(self, mapfn):
#        """ TODO: docstring - mapfn should map old->new variable-index-tuples """
#        new_coeff_dict = { mapfn(k): c for k,c in zip(self.inds,self.coeffs) }
#        self.inds = sorted(list(new_coeff_dict.keys()))
#        self.coeffs = _np.array([new_coeff_dict[k] for k in self.inds], self.coeffs.dtype)
#
#    def mult_poly(self,x):
#        """ Does self * x where x is a polynomial """
#        coeff_dict = {}
#        for k1,v1 in zip(self.inds,self.coeffs):
#            for k2,v2 in zip(x.inds,x.coeffs):
#                k = tuple(sorted(k1+k2))
#                if k in coeff_dict: coeff_dict[k] += v1*v2
#                else: coeff_dict[k] = v1*v2
#        return Polynomial(coeff_dict)
#
#    def mult_scalar(self, x):
#        # assume a scalar that can multiply values
#        newpoly = self.copy()
#        newpoly.coeffs *= x
#        return newpoly
#
#    #def multin_scalar(self, x):
#    #    """ self *= scalar """
#    #    # assume a scalar that can multiply values
#    #    self.coeffs *= x
#
#    def __str__(self):
#        def fmt(x):
#            if abs(_np.imag(x)) > 1e-6:
#                if abs(_np.real(x)) > 1e-6: return "(%.3f+%.3fj)" % (x.real, x.imag)
#                else: return "(%.3fj)" % x.imag
#            else: return "%.3f" % x.real
#            
#        termstrs = []
#        for ik,k in enumerate(self.inds):
#            varstr = ""; last_i = None; n=0
#            for i in sorted(k):
#                if i == last_i: n += 1
#                elif last_i is not None:
#                    varstr += "x%d%s" % (last_i, ("^%d" % n) if n > 1 else "")
#                last_i = i
#            if last_i is not None:
#                varstr += "x%d%s" % (last_i, ("^%d" % n) if n > 1 else "")
#            #print("DB: k = ",k, " varstr = ",varstr)
#            if abs(self.coeffs[ik]) > 1e-4:
#                termstrs.append( "%s%s" % (fmt(self.coeffs[ik]), varstr) )
#        return " + ".join(termstrs)
#
#    def __repr__(self):
#        return "Poly[ " + str(self) + " ]"
#
#    def __add__(self,x):
#        newinds = []
#        newcoeffs = [] # a list for now
#        
#        i1 = i2 = 0 # pointers to self & x, respectively
#        it1 = iter(self.inds)
#        it2 = iter(x.inds)
#
#        try:
#            vinds1 = next(it1)
#            vinds2 = next(it2)
#            while True:
#                if vinds1 == vinds2:
#                    s = self.coeffs[i1] + x.coeffs[i2]
#                    if abs(s) > 1e-12:
#                        newinds.append(vinds1)
#                        newcoeffs.append(s)
#                    i1 += 1; i2 += 1
#                    vinds1 = next(it1)
#                    vinds2 = next(it2)
#                elif vinds1 < vinds2:
#                    newinds.append(vinds1)
#                    newcoeffs.append(self.coeffs[i1])
#                    i1 += 1; vinds1 = next(it1)
#                else:
#                    newinds.append(vinds2)
#                    newcoeffs.append(x.coeffs[i2])
#                    i2 += 1; vinds2 = next(it2)
#        except StopIteration: pass
#        
#        if i1 < len(self.inds):
#            newinds.extend(self.inds[i1:])
#            newcoeffs.extend(self.coeffs[i1:])
#        if i2 < len(x.inds):
#            newinds.extend(x.inds[i2:])
#            newcoeffs.extend(x.coeffs[i2:])
#
#        newpoly = Polynomial()
#        newpoly.inds = newinds
#        newpoly.coeffs = _np.array(newcoeffs, self.coeffs.dtype)
#        return newpoly
#                
#
#    def __mul__(self,x):
#        if isinstance(x, Polynomial):
#            return self.mult_poly(x)
#        else: # assume a scalar that can multiply values
#            return self.mult_scalar(x)
#
#    def __rmul__(self, x):
#        return self.__mul__(x)
#
#    def __pow__(self,n):
#        ret = Polynomial({(): 1.0}) 
#        cur = self
#        for i in range(int(np.floor(np.log2(n)))+1):
#            rem = n % 2 #gets least significant bit (i-th) of n
#            if rem == 1: ret *= cur # add current power of x (2^i) if needed  
#            cur = cur*cur # current power *= 2
#            n //= 2 # shift bits of n right 
#        return ret
#
#    def __copy__(self):
#        return self.copy()


def bulk_eval_compact_polys(compact_poly_tapes, paramvec, dest_shape):
    vtape, ctape = compact_poly_tapes
    result = _np.empty(dest_shape,ctape.dtype) # auto-determine type?
    res = result.flat # for 1D access
    
    c = 0; i = 0; r = 0
    while i < vtape.size:
        poly_val = 0
        nTerms = vtape[i]; i+=1
        #print("POLY w/%d terms (i=%d)" % (nTerms,i))
        for m in range(nTerms):
            nVars = vtape[i]; i+=1 # number of variable indices in this term
            a = ctape[c]; c+=1
            #print("  TERM%d: %d vars, coeff=%s" % (m,nVars,str(a)))
            for k in range(nVars):
                a *= paramvec[ vtape[i] ]; i+=1
            poly_val += a
            #print("  -> added %s to poly_val = %s" % (str(a),str(poly_val))," i=%d, vsize=%d" % (i,vtape.size))
        res[r] = poly_val; r+=1
    assert(c == ctape.size),"Coeff Tape length error: %d != %d !" % (c,ctape.size)
    assert(r == result.size),"Result/Tape size mismatch: only %d result entries filled!" % r
    return result
        
