# Motif.pl
> Discriminatively search for a motif of a given length using given nucleotide sequences
> <br>
> – At least that's what it says on the box...

This is part of a project I worked on in highschool. I thought I'd just throw this here as backup.

The idea was to write a relatively fast and simple [motif finder](http://en.wikipedia.org/wiki/Sequence_motif). The result is a working, semi-accurate command line program that can find motifs in genetic data if you mess with its options for long enough.

## Getting it to work
I don't really remember what packages are needed, but I developed it on Linux a few years ago, and got it working on OS X not so long ago, so given the right packages, this should work on either system.

It seems to require the following packages (though a few more might be needed):

1. `GD::Graph`
1. `Math::GSL`
1. `Text::LevenshteinXS`
 
## The paper
The related paper should be somewhere on the Internet, I'll attach a link when I find it.