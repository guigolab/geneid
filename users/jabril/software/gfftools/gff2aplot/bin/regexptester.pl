#line 1187 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
BEGIN { ($NULL,$T,$F) = ("+@+NULL+@+",1,0); }
    chomp;
    s{^:::\s+}{}o;
    print ">>> $_ \n";
    ($a,$b) = split /\s+:::\s+/o, $_, 2;
    @results = &find_regexp($a);
    @results_mod = (($results[0] ? "TRUE" : "FALSE"),
                    ($results[1] ? "TRUE" : "FALSE"),
                    (defined($results[2]) ? $results[2] : "FALSE"),
                    "$results[3]");
    printf ">>> IS_OK::%s::  IS_NEG::%s::  GET_ID::%s::  GET_STRING::%s::\n",
           @results_mod;
    !$results[1] &&
         print ">>> ".
             (  ( (eval { $b =~ m{$results[3]}; 1 } || 0) ?
                         ($b =~ m{$results[3]}) : 0 ) ? 
                     "\"$b\" MATCHES \/$results[3]\/ " : 
                     "\"$b\" DOES NOT MATCH \/$results[3]\/ ")."\n\n";
    $results[1] &&
         print ">>> ".
             (!(( (eval { $b =~ m{$results[3]}; 1 } || 0) ?
                         ($b =~ m{$results[3]}) : 0 )) ? 
                       "\"$b\"  MATCHES !\/$results[3]\/ " : 
                       "\"$b\" DOES NOT MATCH !\/$results[3]\/ ")."\n\n";
sub find_regexp() {
    my $string = $_[0];
    my ($isOK_flg,$not_flg,$id_flg,$tmpstr,$tmpid);
    $isOK_flg = $T;
    $not_flg = $F;
    $id_flg = undef;
    $string =~ s{^!}{}o && ($not_flg = $T); # not_regexp is true
    $string =~ s{(\\@)$}{@@}o;
    $string = &escape_input($string);
    ($tmpstr, $tmpid) = (undef, undef);
    ( reverse($string) =~ m{^([^\/@]*?)(?:@){1}(.*)$}o ) && do {
        $tmpstr = reverse($2);
        $tmpid  = reverse($1);
    };
  REGEXPS: {
      print ">>> STRING($tmpstr) ID($tmpid)\n";
      (defined($tmpid) && $tmpid ne "") && ($id_flg = $tmpid);
      (defined($tmpstr) && $tmpstr ne "") || do {
          $string eq '@' && ($string="", $isOK_flg=$F);
          $tmpstr = $string;
      };
      $tmpstr eq '*' && do {
          $string = '^.*$';
          last REGEXPS;
      };
      $tmpstr =~ m{^/(.*)/$}o && do {    
          ($string, $isOK_flg) = &eval_regexp($1);
          last REGEXPS;
      }; # $tmpstr is a regexp
      $string = '^'.(quotemeta($tmpstr)).'$';
    }; # REGEXPS
    return ($isOK_flg, $not_flg, $id_flg, $string);
} # find_regexp
sub eval_regexp() {
    my $str = $_[0];
    my $flag;
    # $str =~ s{^[^\\]*?/+}{}o;
    # $str =~ s{/+.*?$}{}o;
    eval { "" =~ m{$str}; $flag = $T; } || ($flag = $F);
    return ($str, $flag);
} # eval_regexp
sub escape_input() {
  $_[0] =~ s{([;<>&!{}'`"])}{\\$1}og;
  return $_[0];
} # escape_input
