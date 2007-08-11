package gphase.regregex;
public class RegTest
{
    /*
     * Copyright (c) 2005, Damien Mascord <tusker@tusker.org> All rights reserved.
     * 
     * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following
     * conditions are met:
     * 
     * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
     * disclaimer in the documentation and/or other materials provided with the distribution. Neither the name of the <ORGANIZATION>
     * nor the names of its contributors may be used to endorse or promote products derived from this software without specific
     * prior written permission. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
     * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
     * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
     * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
     * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
     * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
     * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
     * 
     */

    private static final String[] _re = {"^(([^:]+)://)?([^:/]+)(:([0-9]+))?(/.*)", // URL match
        "(([^:]+)://)?([^:/]+)(:([0-9]+))?(/.*)", // URL match without starting ^
        "usd [+-]?[0-9]+.[0-9][0-9]", // Canonical US dollar amount
        "\\b(\\w+)(\\s+\\1)+\\b", // Duplicate words
        "\\{(\\d+):(([^}](?!-} ))*)", // this is meant to match against the "some more text and ..." but it causes ORO Matcher
    // to fail, so we won't include this by default... it is also WAY too slow to test
    // we will test [4][5] only 1 iteration...
    };

    private static final String[] _str = {
        "http://www.linux.com/",
        "http://www.thelinuxshow.com/main.php3",
        "usd 1234.00",
        "he said she said he said no",
        "same same same",
        "{1:\n" + "this is some more text - and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more\n"
            + "this is some more text and some more and some more and even more at the end\n" + "-}\n", // very large bit of text...

    };

    private static boolean[][] expectedMatch = new boolean[_re.length][_str.length];

    static
    {
        expectedMatch[0][0] = true;
        expectedMatch[0][1] = true;
        expectedMatch[0][2] = false;
        expectedMatch[0][3] = false;
        expectedMatch[0][4] = false;
        expectedMatch[0][5] = false;
        expectedMatch[1][0] = true;
        expectedMatch[1][1] = true;
        expectedMatch[1][2] = false;
        expectedMatch[1][3] = false;
        expectedMatch[1][4] = false;
        expectedMatch[1][5] = false;
        expectedMatch[2][0] = false;
        expectedMatch[2][1] = false;
        expectedMatch[2][2] = true;
        expectedMatch[2][3] = false;
        expectedMatch[2][4] = false;
        expectedMatch[2][5] = false;
        expectedMatch[3][0] = false;
        expectedMatch[3][1] = false;
        expectedMatch[3][2] = false;
        expectedMatch[3][3] = false;
        expectedMatch[3][4] = true;
        expectedMatch[3][5] = false;
        expectedMatch[4][0] = false;
        expectedMatch[4][1] = false;
        expectedMatch[4][2] = false;
        expectedMatch[4][3] = false;
        expectedMatch[4][4] = true;
        expectedMatch[4][5] = false;
    }

    private static boolean debug = true;
    private static boolean html = false;

    private final static int ITERATIONS = 10000;

    public static final void main(String[] args)
    {
	    try
	    {
	        // org.apache.regexp.* test
	        long[][][] timeTaken = new long[_re.length][ITERATIONS][_str.length];
	        boolean[][] matches = new boolean[_re.length][_str.length];
	        long startTime = System.currentTimeMillis();
	        // ----------------------//

            // RegularExpression.RE version [fails on URL tests]
            /*
             * startTime = System.currentTimeMillis(); for (int regnum = 0; regnum < _re.length; regnum++) { RegularExpression.RE
             * regexpr = new RegularExpression.RE(_re[regnum], true);
             * 
             * for (int itter=0; itter<ITERATIONS; itter++){ long iterStarTime = System.currentTimeMillis(); for (int
             * strnum=0;strnum<_str.length;strnum++){ boolean b = regexpr.matches(_str[strnum]); timeTaken[regnum][itter][strnum] = (
             * System.currentTimeMillis() - iterStarTime ) ; if (debug && itter==0){ System.out.println(_re[regnum] + " against " +
             * _str[strnum] + ":" + b); } } } } endTime = System.currentTimeMillis(); printResult("RegularExpression.RE", timeTaken,
             * (endTime - startTime), html );
             */

            // ----------------------//
            // gnu.rex.Rex version [too much crap is printed to screen]
            /*
             * startTime = System.currentTimeMillis(); for (int regnum = 0; regnum < _re.length; regnum++) {
             * gnu.rex.Rex.config_GroupBraces("(", ")"); gnu.rex.Rex.config_Alternative("|"); gnu.rex.Rex regexpr =
             * gnu.rex.Rex.build(_re[regnum]); for (int itter=0; itter<ITERATIONS; itter++){ long iterStarTime =
             * System.currentTimeMillis(); for (int strnum=0;strnum<_str.length;strnum++){ gnu.rex.RexResult rexResult =
             * regexpr.match(_str[strnum].toCharArray(),0,0); timeTaken[regnum][itter][strnum] = ( System.currentTimeMillis() -
             * iterStarTime ) ; } } } endTime = System.currentTimeMillis(); System.out.println("gnu.rex.Rex took " + (endTime -
             * startTime) + "ms");
             */
            // ----------------------//
            // dk.brics.automaton.RegExp version [fails on URL tests]
            startTime = System.currentTimeMillis();
            for (int regnum = 0; regnum < _re.length; regnum++)
            {
                try
                {
                    RegExp regexpr = new RegExp(_re[regnum]);
                    Automaton auto = regexpr.toAutomaton();
                    RunAutomaton runauto = new RunAutomaton(auto, true);
                    for (int itter = 0; itter < ITERATIONS; itter++)
                    {
                        long iterStarTime = System.currentTimeMillis();
                        for (int strnum = 0; strnum < _str.length; strnum++)
                        {
                            boolean b = runauto.run(_str[strnum]);
                            matches[regnum][strnum] = (b == expectedMatch[regnum][strnum]);
                            timeTaken[regnum][itter][strnum] = (System.currentTimeMillis() - iterStarTime);
                            if (debug && itter == 0)
                            {
                                System.out.println(_re[regnum] + " against " + _str[strnum] + ":" + b);
                            }
                            // only test the big one 1 iteration only
                            if (strnum == 5)
                            {
                                break;
                            }
                        }
                        // only test the big one 1 iteration only
                        if (regnum == 4)
                        {
                            break;
                        }
                    }
                }
                catch (Throwable e)
                {
                    if (debug)
                    {
                        System.out.println(_re[regnum] + "  failed badly");
                    }
                }
            }
            long endTime = System.currentTimeMillis();
            printResult("dk.brics.automaton.RegExp", timeTaken, (endTime - startTime), matches, html);
        } catch (Exception e)
        {
            e.printStackTrace();
        }
    }

	private static final void printResult(String regexName, long[][][] matrix, long totalTime, boolean[][] matches, boolean html)
    {
        // timeTaken[regnum][itter][strnum]
        if (html)
        {
            System.out.println("<table>");
            System.out.println("<tr><th colspan=\"3\"><h2>Regular expression library:</h2></th><td colspan=\"3\"><h2>" + regexName
                + "</h2></td></tr>");
        }
        else
        {
            System.out.println("------------------------------------------");
            System.out.println("Regular expression library: " + regexName + "\n");
        }
        for (int re = 0; re < _re.length; re++)
        {
            if (html)
            {
                System.out.println("<tr><th>RE:</th><td colspan=\"5\">" + _re[re] + "</td></tr>");
                System.out
                    .println("<tr><th>MS</th><th>MAX</th><th>AVG</th><th>MIN</th><th>DEV</th><th>INPUT</th><th>MATCH</th></tr>");
            }
            else
            {
                System.out.println("RE: " + _re[re]);
                System.out.println("  MS\tMAX\tAVG\tMIN\tDEV\tINPUT\tMATCH");
            }
            for (int str = 0; str < _str.length; str++)
            {
                long total = 0;
                long sumOfSq = 0;
                long min = Long.MAX_VALUE;
                long max = Long.MIN_VALUE;
                for (int i = 0; i < ITERATIONS; i++)
                {
                    long elapsed = matrix[re][i][str];
                    total += elapsed;
                    sumOfSq += elapsed * elapsed;
                    if (elapsed < min)
                    {
                        min = elapsed;
                    }
                    if (elapsed > max)
                    {
                        max = elapsed;
                    }
                }
                // calc std dev
                long stdDev = (long) java.lang.Math.sqrt((sumOfSq - ((total * total) / ITERATIONS)) / (ITERATIONS - 1));

                if (html)
                {
                    System.out.println("<tr><td>" + total + "</td><td>" + max + "</td><td>" + (double) total / ITERATIONS
                        + "</td><td>" + min + "</td><td>" + stdDev + "</td><td>" + _str[str] + "</td><td>" + matches[re][str]
                        + "</td></tr>");
                }
                else
                {
                    System.out.println("  " + total + "\t" + max + "\t" + (double) total / ITERATIONS + "\t" + min + "\t" + stdDev
                        + "\t'" + _str[str] + "\t'" + matches[re][str] + "'");
                }
            }
        }
        if (html)
        {
            System.out.println("<tr><th colspan=\"3\"><h2>Total time taken:</h2></th><td colspan=\"3\"><h2>" + totalTime
                + "</h2></td></tr>");
            System.out.println("</table>");
        }
        else
        {
            System.out.println("Total time taken: " + totalTime);
            System.out.println("------------------------------------------");
        }
    }

}
