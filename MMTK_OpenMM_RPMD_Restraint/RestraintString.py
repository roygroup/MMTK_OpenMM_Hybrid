#This function creates an energy expression as a string based on 2 inputted list of atoms, 
#spring constant, and equilibrium value as well as what type of expression as a string,
#currently only 'com' for centre of mass and 'cog' for centre of geometry
#Last modified: June 17, 2013

"""
Copyright (c) 2013 Stephen Constable, Nabil Faruk, Kevin Bishop, Pierre-Nicholas Roy
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

def CreateRestraintString(g1,g2,k0,r0,option):

    m = 'm'
    x = 'x'
    y = 'y'
    z = 'z'

    if option == 'cog':
        print 'Creating Centre of Geometry Restraint Energy Expression:'

        x1=y1=z1=M1 = '('
        x2=y2=z2=M2 = '('

        for i in range(0,len(g1)):
            if i != 0:
                x1 += ' + '
                y1 += ' + '
                z1 += ' + '
            k = str(i+1)
            x1 += x + k
            y1 += y + k
            z1 += z + k

        x1 += ')'
        y1 += ')'
        z1 += ')'
        M1 = str(len(g1))

        for i in range(len(g1),len(g1)+len(g2)):
            if i != len(g1):
                x2 += ' + '
                y2 += ' + '
                z2 += ' + '
            k = str(i+1)
            x2 += x + k
            y2 += y + k
            z2 += z + k

        x2 += ')'
        y2 += ')'
        z2 += ')'
        M2 = str(len(g2))

        dx = '((' + x1 + ' / ' + M1 + ') - (' + x2 + ' / ' + M2 +  '))' 
        dy = '((' + y1 + ' / ' + M1 + ') - (' + y2 + ' / ' + M2 +  '))' 
        dz = '((' + z1 + ' / ' + M1 + ') - (' + z2 + ' / ' + M2 +  '))'

        dist = 'sqrt(' + dx + '^2 + ' + dy + '^2 + ' + dz + '^2)'

        return_string = '0.5 * ' +  str(float(k0)) + ' * '+ '(' + str(float(r0)) + ' - ' + dist + ')^2'

        return return_string

    elif option == 'com':
        print 'Creating Centre of Mass Restraint Energy Expression:'

        x1=y1=z1=M1= '('
        x2=y2=z2=M2= '('

        for i in range(0,len(g1)):
            if i != 0:
                x1 += ' + '
                y1 += ' + '
                z1 += ' + '
                M1 += ' + '
            k = str(i+1)
            x1 += m + k + '*' + x + k
            y1 += m + k + '*' + y + k
            z1 += m + k + '*' + z + k
            M1 += m + k

        x1 += ')'
        y1 += ')'
        z1 += ')'
        M1 += ')'

        for i in range(len(g1),len(g1)+len(g2)):
            if i != len(g1):
                x2 += ' + '
                y2 += ' + '
                z2 += ' + '
                M2 += ' + '
            k = str(i+1)
            x2 += m + k + '*' + x + k
            y2 += m + k + '*' + y + k
            z2 += m + k + '*' + z + k
            M2 += m + k

        x2 += ')'
        y2 += ')'
        z2 += ')'
        M2 += ')'

        dx = '((' + x1 + ' / ' + M1 + ') - (' + x2 + ' / ' + M2 +  '))' 
        dy = '((' + y1 + ' / ' + M1 + ') - (' + y2 + ' / ' + M2 +  '))' 
        dz = '((' + z1 + ' / ' + M1 + ') - (' + z2 + ' / ' + M2 +  '))' 

        dist = 'sqrt(' + dx + '^2 + ' + dy + '^2 + ' + dz + '^2)'

        return_string = '0.5 * ' +  str(float(k0)) + ' * '+ '(' + str(float(r0)) + ' - ' + dist + ')^2'

        return return_string

    else:
        print 'Unrecognized Restraint Option'

