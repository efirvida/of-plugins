#ifndef UTILITIES_H
#define UTILITIES_H

// Filter-out debug messages when not required
#ifdef ADAPTER_DEBUG_MODE
#define DEBUG(x) x
#else
#define DEBUG(x)
#endif

// Filter-out timing when not required
#ifdef ADAPTER_ENABLE_TIMINGS
#define TIMING_MODE(x) x
#else
#define TIMING_MODE(x)
#endif

// Shortcut for constructing clocks
#ifdef ADAPTER_ENABLE_TIMINGS
#define SETUP_TIMER() \
    clockValue clock; \
    clock.update()
#define REUSE_TIMER() clock.update()
#else
#define SETUP_TIMER()
#define REUSE_TIMER()
#endif


// Shortcut for adding time
#ifdef ADAPTER_ENABLE_TIMINGS
#define ACCUMULATE_TIMER(x) x += clock.elapsed()
#else
#define ACCUMULATE_TIMER(x)
#endif

// String added in the beginning of every printed message
#define INFO_STR_ADAPTER "---[preciceAdapter] "

#include "IOstreams.H"
#include "Pstream.H"
#include "IPstream.H"
#include "OPstream.H"

namespace preciceAdapter
{

//- Gathers a single value from all ranks into a List on the master rank.
template<class T>
void gather(Foam::List<T>& allRanksValues)
{
    if (!Foam::Pstream::parRun()) return;

    if (!Foam::Pstream::master())
    {
        Foam::OPstream toMaster(Foam::Pstream::commsTypes::scheduled, Foam::Pstream::masterNo());
        toMaster << allRanksValues[Foam::Pstream::myProcNo()];
    }
    else
    {
        for (int p = 1; p < Foam::Pstream::nProcs(); ++p)
        {
            Foam::IPstream fromProc(Foam::Pstream::commsTypes::scheduled, p);
            fromProc >> allRanksValues[p];
        }
    }
}

//- Broadcasts a value from the master rank to all other ranks.
template<class T>
void broadcast(T& value)
{
    if (!Foam::Pstream::parRun()) return;
    Foam::Pstream::broadcast(value);
}

//- Scatters a List of values from the master rank into individual values on each rank.
template<class T>
void scatter(Foam::List<T>& allRanksValues)
{
    if (!Foam::Pstream::parRun()) return;

    if (Foam::Pstream::master())
    {
        for (int p = 1; p < Foam::Pstream::nProcs(); ++p)
        {
            Foam::OPstream toProc(Foam::Pstream::commsTypes::scheduled, p);
            toProc << allRanksValues[p];
        }
    }
    else
    {
        Foam::IPstream fromMaster(Foam::Pstream::commsTypes::scheduled, Foam::Pstream::masterNo());
        fromMaster >> allRanksValues[Foam::Pstream::myProcNo()];
    }
}

} // End namespace preciceAdapter

void adapterInfo(const std::string message, const std::string level = "debug");

// Returns true if the match is found in the beginning of the string s, ignoring case
bool matchingStrings(std::string s, std::string match);

#endif
