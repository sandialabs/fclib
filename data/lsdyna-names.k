$ Copyright (2000) Sandia Corporation. Under the terms of Contract
$ DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains 
$ certain rights in this software.
$ 
$ Redistribution and use in source and binary forms, with or without
$ modification, are permitted provided that the following conditions are met:
$
$   * Redistributions of source code must retain the above copyright
$     notice, this list of conditions and the following disclaimer.
$
$   * Redistributions in binary form must reproduce the above copyright
$     notice, this list of conditions and the following disclaimer in the
$     documentation and/or other materials provided with the
$     distribution.
$
$   * Neither the name of Sandia nor the names of any contributors may
$     be used to endorse or promote products derived from this software
$     without specific prior written permission.
$
$   * Modified source versions must be plainly marked as such, and must
$     not be misrepresented as being the original software.
$
$ THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
$ AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
$ IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
$ ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR 
$ ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
$ DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
$ SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
$ CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
$ LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
$ OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH 
$ DAMAGE.

$ This is not a real LSDya input file.
$ This was fabricated to test various input file parsers.
$ $Date$
*KEYWORD

*PARAMETER
I idfor21 21
*PARAMETER
$ could have a comment here
R idfor22 22
*PARAMETER
$ could have
$ multiple
$ comments
I idfor23 23
*PARAMETER
$ saw one once like this
Iidfor24 24
*PARAMETER
I idfor25 25
*PARAMETER
I idfor26 26

$ easy part declaration
*PART
part1
1,1,1,1
*PART
$    PID   SECID  MID   <- comment before
part2
       2       1    1
*PART
$    PID   SECID  MID   <- comments before
$    PID   SECID  MID   <- comments before
part3
       3       1    1
*PART
part4
$    PID   SECID  MID   <- comment between
       4       1    1
*PART
part5
$    PID   SECID  MID   <- comments between
$    PID   SECID  MID   <- comments between
       5       1    1
*PART
$    PID   SECID  MID   <- comments all over
part6
$    PID   SECID  MID   <- comments all over
      6       1    1

$ part dedlarations w/ parameters$ easy part declaration
*PART
part21
&idfor21,1,1,1
*PART
$    PID   SECID  MID   <- comment before
part22
 &idfor22       1    1
*PART
$    PID   SECID  MID   <- comments before
$    PID   SECID  MID   <- comments before
part23
 &idfor23       1    1
*PART
part24
$    PID   SECID  MID   <- comment between
 &idfor24       1    1
*PART
part25
$    PID   SECID  MID   <- comments between
$    PID   SECID  MID   <- comments between
 &idfor25       1    1
*PART
$    PID   SECID  MID   <- comments all over
part26
$    PID   SECID  MID   <- comments all over
 &idfor26       1    1

$ Try different kinds of names
*PART
part-31-ohno
31,1,1,1
*PART
part_32_ohyes
32,1,1,1
*PART
part.33.me.oh.my
33,1,1,1
*PART
part 34 where are you?
34,1,1,1
$ spaces at front
*PART
  part 35
35,1,1,1
$ spaces at end
*PART
part 36  
36,1,1,1

$ other "bugs"
$ This should be treated as a part
*PART_contact
part41
41,1,1,1
$ Had a bug - single character names didn't work
*PART
a
42,1,1,1
