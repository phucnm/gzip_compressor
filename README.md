# gzip_compressor
A bs gzip compressor

This's a simple gzip compressor compliant to standard `gzip` tool on Linux created for educational purpose. In other words, `gzip` can perfectly unzip compressed files generated from this compressor without warnings. This program accepts input from standard input and produce compressed bitstream to standard output. Implementation was followed [`gzip`](https://tools.ietf.org/html/rfc1952) and [`DEFLATE` ](https://tools.ietf.org/html/rfc1951) specification.

## Usage 

- Build by running `make`
- Run the program by running `./uvgz < inputfile > inputfile.gz` where `inputfile` is an input file we want to compress regardless of text file or binary file.
- To decompress with `gzip` running `gzip -d < inputfile.gz > inputfile`

## Some implementation notes
- This compressor use the maximum window size allowed by `gzip`: 32768 bytes.
- Block header for a block of type 2 is efficiently encoded by a RLE step as described in section 3.2.7. of RFC-1951.
- A simple back-reference algorithm was implemented using a hash table. Essentially, there's a hash chain to keep reference to every group of 3 chars. Thus a reference lookup only costs O(1) time. In the worst case where a window contains only 1 character e.g. "aaaaaa", it can cost O(n) time. A maximum number of chain hit of 512 was used to mitigate these kinds of inputs. This setting should be dynamically configured as an compression parameter. A smaller number can increase compression speed but the compression ratio could be less optimal.
- Package-Merge algorithm was implemented to generate code lengths for literal/length codes, distance codes and code length codes instead of Huffman coding. The purpose is to limit the maximum number of code length value and to generate a more optimal set of lengths from a set of frequencies. The `package_merge` method inspired from an example at page 54 in the book [Intro to Data Compression](https://www.sciencedirect.com/book/9780126208627/introduction-to-data-compression).

