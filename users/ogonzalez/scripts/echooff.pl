Password How to turn echo off for stdin in Perl 
                        


                      use Term::ReadKey;

                      print "Enter your password: ";
                      ReadMode 'noecho';
                      $password = ReadLine 0;
                      ReadMode 'normal';
