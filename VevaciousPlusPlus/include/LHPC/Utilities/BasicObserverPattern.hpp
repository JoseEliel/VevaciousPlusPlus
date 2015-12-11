/*
 * BasicObserverPattern.hpp
 *
 *  Created on: Nov 26, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 *
 *      This contains the two classes needed for a very basic observer pattern:
 *      BasicObserver and BasicObserved. Classes which inherit from
 *      BasicObserver need only to implement RespondToObservedSignal(), while
 *      classes which inherit from BasicObserved don't need to implement
 *      anything.
 */

#ifndef BASICOBSERVERPATTERN_HPP_
#define BASICOBSERVERPATTERN_HPP_

#include <list>

namespace LHPC
{
  // These are the two classes needed for a very basic observer pattern:
  // BasicObserver and BasicObserved. Classes which inherit from BasicObserver
  // need only to implement RespondToObservedSignal(), while classes which
  // inherit from BasicObserved don't need to implement anything.

  // This is an abstract base class that allows BasicObserved objects to call
  // RespondToObservedSignal() on its observers.
  class BasicObserver
  {
  public:
    // The constructor just initializes the list.
    BasicObserver() : stillObservingFlags() {}

    // The destructor needs to inform any observed targets still around to stop
    // to signaling this observer.
    virtual ~BasicObserver() { ClearObservingFlags(); }

    // This should be the function to perform the work which should happen on
    // receiving a signal from the observed object.
    virtual void RespondToObservedSignal() = 0;

    void AcceptFlagFromObserved( bool* const flagFromObserved )
    { stillObservingFlags.push_back( flagFromObserved ); }

    void DiscardFlagFromObserved( bool* const flagFromObserved )
    { stillObservingFlags.remove( flagFromObserved ); }

    void ClearObservingFlags();


  protected:
    std::list< bool* > stillObservingFlags;
  };

  // This class holds a list of BasicObserver pointers and calls
  // RespondToObservedSignal() on them all with UpdateObservers().
  class BasicObserved
  {
  public:
    BasicObserved() : observerList(),
                      observerIterator() {}
    virtual ~BasicObserved() { RemoveAllObservers(); }

    // This goes through observerList and calls its RespondToObservedSignal()
    // if its bool is true, otherwise it removes the observer.
    void UpdateObservers();

    // This registers joiningObserver as an observer which will have its
    // RespondToObservedSignal() function called when this BasicObserved calls
    // its UpdateObservers().
    void RegisterObserver( BasicObserver* const joiningObserver );

    // This removes leavingObserver from the list of observers, also removing
    // any observers in the list which flagged their bools as false, as would
    // happen if their destructors were called.
    void RemoveObserver( BasicObserver* const leavingObserver );

    // This clears observerList after asking all its observers to discard their
    // bool pointers for this BasicObserved instance.
    void RemoveAllObservers();


  protected:
    std::list< std::pair< BasicObserver*, bool > > observerList;
    std::list< std::pair< BasicObserver*, bool > >::iterator observerIterator;
  };





  inline void BasicObserver::ClearObservingFlags()
  {
    for( std::list< bool* >::iterator
         flagIterator( stillObservingFlags.begin() );
         stillObservingFlags.end() != flagIterator;
         ++flagIterator )
    {
      *(*flagIterator) = false;
    }
    stillObservingFlags.clear();
  }





  // This goes through observerList and calls its RespondToObservedSignal()
  // if its bool is true, otherwise it removes the observer.
  inline void BasicObserved::UpdateObservers()
  {
    observerIterator = observerList.begin();
    while( observerList.end() != observerIterator )
    {
      if( observerIterator->second )
      {
        observerIterator->first->RespondToObservedSignal();
        ++observerIterator;
      }
      else
      {
        observerIterator = observerList.erase( observerIterator );
      }
    }
  }

  // This registers joiningObserver as an observer which will have its
  // RespondToObservedSignal() function called when this BasicObserved calls
  // its UpdateObservers().
  inline void
  BasicObserved::RegisterObserver( BasicObserver* const joiningObserver )
  {
    observerList.push_back( std::pair< BasicObserver*, bool >( joiningObserver,
                                                               true ) );

    // New observers are given a pointer to their associated bools, which
    // indicate that their observers should be fine for
    // RespondToObservedSignal() calls if the bool is true, but if the bool is
    // false, this BasicObserved instance knows to remove the observer with its
    // next pass through the list.
    joiningObserver->AcceptFlagFromObserved( &(observerList.back().second) );
  }

  // This removes leavingObserver from the list of observers, also removing any
  // observers in the list which flagged their bools as false, as would happen
  // if their destructors were called.
  inline void
  BasicObserved::RemoveObserver( BasicObserver* const leavingObserver )
  {
    observerIterator = observerList.begin();
    while( observerList.end() != observerIterator )
    {
      if( !(observerIterator->second) )
      {
        observerIterator = observerList.erase( observerIterator );
      }
      else if( leavingObserver == observerIterator->first )
      {
        leavingObserver->DiscardFlagFromObserved(
                                                 &(observerIterator->second) );
        observerIterator = observerList.erase( observerIterator );
      }
      else
      {
        ++observerIterator;
      }
    }
  }

  // This clears observerList after asking all its observers to discard their
  // bool pointers for this BasicObserved instance.
  inline void BasicObserved::RemoveAllObservers()
  {
    observerIterator = observerList.begin();
    while( observerList.end() != observerIterator )
    {
      if( observerIterator->second )
      {
        observerIterator->first->DiscardFlagFromObserved(
                                                 &(observerIterator->second) );
      }
      ++observerIterator;
    }
    observerList.clear();
  }

}

#endif /* BASICOBSERVERPATTERN_HPP_ */
