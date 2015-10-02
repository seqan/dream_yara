//![include]
#include <seqan/bam_io.h>

using namespace seqan;

int main()
{
//![include]
//![ctor]
    // Open input BAM file, BamFileIn supports both SAM and BAM files.
    BamFileIn bamFileIn("example.bam");
//![ctor]

//![open]
    // Open output SAM file by passing the context of bamFileIn and the filename to open.
    BamFileOut samFileOut(context(bamFileIn), toCString(samFileOutName));
//![open]

//![header]
    // Copy header.
    BamHeader header;
    readHeader(header, bamFileIn);
    writeHeader(samFileOut, header);
//![header]

//![records]
    // Copy all records.
    BamAlignmentRecord record;
    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);
        writeRecord(samFileOut, record);
    }

    return 0;
}
//![records]
