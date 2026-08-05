# 1 "mpp/mpp.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "mpp/mpp.F90"
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!-----------------------------------------------------------------------
!                 Communication for message-passing codes
!
! AUTHOR: V. Balaji (V.Balaji@noaa.gov)
!         SGI/GFDL Princeton University
!
!-----------------------------------------------------------------------

!> @defgroup mpp_mod mpp_mod
!! @ingroup mpp
!! @{
!! @brief This module defines interfaces for common operations using message-passing libraries.
!! Any type-less arguments in the documentation are MPP_TYPE_ which is defined by the pre-processor
!! to create multiple subroutines out of one implementation for use in an interface. See the note
!! below for more information
!!
!! @author V. Balaji <"V.Balaji@noaa.gov">
!!
!!   A set of simple calls to provide a uniform interface
!!   to different message-passing libraries. It currently can be
!!   implemented either in the SGI/Cray native SHMEM library or in the MPI
!!   standard. Other libraries (e.g MPI-2, Co-Array Fortran) can be
!!   incorporated as the need arises.
!!
!!   The data transfer between a processor and its own memory is based
!!   on <TT>load</TT> and <TT>store</TT> operations upon
!!   memory. Shared-memory systems (including distributed shared memory
!!   systems) have a single address space and any processor can acquire any
!!   data within the memory by <TT>load</TT> and
!!   <TT>store</TT>. The situation is different for distributed
!!   parallel systems. Specialized MPP systems such as the T3E can simulate
!!   shared-memory by direct data acquisition from remote memory. But if
!!   the parallel code is distributed across a cluster, or across the Net,
!!   messages must be sent and received using the protocols for
!!   long-distance communication, such as TCP/IP. This requires a
!!   ``handshaking`` between nodes of the distributed system. One can think
!!   of the two different methods as involving <TT>put</TT>s or
!!   <TT>get</TT>s (e.g the SHMEM library), or in the case of
!!   negotiated communication (e.g MPI), <TT>send</TT>s and
!!   <TT>recv</TT>s.
!!
!!   The difference between SHMEM and MPI is that SHMEM uses one-sided
!!   communication, which can have very low-latency high-bandwidth
!!   implementations on tightly coupled systems. MPI is a standard
!!   developed for distributed computing across loosely-coupled systems,
!!   and therefore incurs a software penalty for negotiating the
!!   communication. It is however an open industry standard whereas SHMEM
!!   is a proprietary interface. Besides, the <TT>put</TT>s or
!!   <TT>get</TT>s on which it is based cannot currently be implemented in
!!   a cluster environment (there are recent announcements from Compaq that
!!   occasion hope).
!!
!!   The message-passing requirements of climate and weather codes can be
!!   reduced to a fairly simple minimal set, which is easily implemented in
!!   any message-passing API. <TT>mpp_mod</TT> provides this API.
!!
!!    Features of <TT>mpp_mod</TT> include:
!!    <ol>
!!     <li> Simple, minimal API, with free access to underlying API for </li>
!!       more complicated stuff.<BR/>
!!     <li> Design toward typical use in climate/weather CFD codes. </li>
!!     <li> Performance to be not significantly lower than any native API. </li>
!!    </ol>
!!
!!   This module is used to develop higher-level calls for
!!   domain decomposition (@ref mpp_domains) and parallel I/O (@ref fms2_io)
!! <br/>
!!   Parallel computing is initially daunting, but it soon becomes
!!   second nature, much the way many of us can now write vector code
!!   without much effort. The key insight required while reading and
!!   writing parallel code is in arriving at a mental grasp of several
!!   independent parallel execution streams through the same code (the SPMD
!!   model). Each variable you examine may have different values for each
!!   stream, the processor ID being an obvious example. Subroutines and
!!   function calls are particularly subtle, since it is not always obvious
!!   from looking at a call what synchronization between execution streams
!!   it implies. An example of erroneous code would be a global barrier
!!   call (see @ref mpp_sync below) placed
!!   within a code block that not all PEs will execute, e.g:
!!
!!   <PRE>
!!   if( pe.EQ.0 )call mpp_sync()
!!   </PRE>
!!
!!   Here only PE 0 reaches the barrier, where it will wait
!!   indefinitely. While this is a particularly egregious example to
!!   illustrate the coding flaw, more subtle versions of the same are
!!   among the most common errors in parallel code.
!!  <br/>
!!   It is therefore important to be conscious of the context of a
!!   subroutine or function call, and the implied synchronization. There
!!   are certain calls here (e.g <TT>mpp_declare_pelist, mpp_init,
!!   mpp_set_stack_size</TT>) which must be called by all
!!   PEs. There are others which must be called by a subset of PEs (here
!!   called a <TT>pelist</TT>) which must be called by all the PEs in the
!!   <TT>pelist</TT> (e.g <TT>mpp_max, mpp_sum, mpp_sync</TT>). Still
!!   others imply no synchronization at all. I will make every effort to
!!   highlight the context of each call in the MPP modules, so that the
!!   implicit synchronization is spelt out.
!! <br/>
!!   For performance it is necessary to keep synchronization as limited
!!   as the algorithm being implemented will allow. For instance, a single
!!   message between two PEs should only imply synchronization across the
!!   PEs in question. A <I>global</I> synchronization (or <I>barrier</I>)
!!   is likely to be slow, and is best avoided. But codes first
!!   parallelized on a Cray T3E tend to have many global syncs, as very
!!   fast barriers were implemented there in hardware.
!! <br/>
!!   Another reason to use pelists is to run a single program in MPMD
!!   mode, where different PE subsets work on different portions of the
!!   code. A typical example is to assign an ocean model and atmosphere
!!   model to different PE subsets, and couple them concurrently instead of
!!   running them serially. The MPP module provides the notion of a
!!   <I>current pelist</I>, which is set when a group of PEs branch off
!!   into a subset. Subsequent calls that omit the <TT>pelist</TT> optional
!!   argument (seen below in many of the individual calls) assume that the
!!   implied synchronization is across the current pelist. The calls
!!   <TT>mpp_root_pe</TT> and <TT>mpp_npes</TT> also return the values
!!   appropriate to the current pelist. The <TT>mpp_set_current_pelist</TT>
!!   call is provided to set the current pelist.
!! </DESCRIPTION>
!! <br/>
!!
!!  @note F90 is a strictly-typed language, and the syntax pass of the
!!  compiler requires matching of type, kind and rank (TKR). Most calls
!!  listed here use a generic type, shown here as <TT>MPP_TYPE_</TT>. This
!!  is resolved in the pre-processor stage to any of a variety of
!!  types. In general the MPP operations work on 4-byte and 8-byte
!!  variants of <TT>integer, real, complex, logical</TT> variables, of
!!  rank 0 to 5, leading to 48 specific module procedures under the same
!!  generic interface. Any of the variables below shown as
!!  <TT>MPP_TYPE_</TT> is treated in this way.

module mpp_mod

! Define rank(X) for PGI compiler








  use gfdl_nompi_f08


  use iso_fortran_env,   only : INPUT_UNIT, OUTPUT_UNIT, ERROR_UNIT
  use mpp_parameter_mod, only : MPP_VERBOSE, MPP_DEBUG, ALL_PES, ANY_PE, NULL_PE
  use mpp_parameter_mod, only : NOTE, WARNING, FATAL, MPP_CLOCK_DETAILED,MPP_CLOCK_SYNC
  use mpp_parameter_mod, only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE_DRIVER
  use mpp_parameter_mod, only : CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA
  use mpp_parameter_mod, only : MAX_EVENTS, MAX_BINS, MAX_EVENT_TYPES, MAX_CLOCKS
  use mpp_parameter_mod, only : MAXPES, EVENT_WAIT, EVENT_ALLREDUCE, EVENT_BROADCAST
  use mpp_parameter_mod, only : EVENT_ALLTOALL
  use mpp_parameter_mod, only : EVENT_TYPE_CREATE, EVENT_TYPE_FREE
  use mpp_parameter_mod, only : EVENT_RECV, EVENT_SEND, MPP_READY, MPP_WAIT
  use mpp_parameter_mod, only : mpp_parameter_version=>version
  use mpp_parameter_mod, only : DEFAULT_TAG
  use mpp_parameter_mod, only : COMM_TAG_1,  COMM_TAG_2,  COMM_TAG_3,  COMM_TAG_4
  use mpp_parameter_mod, only : COMM_TAG_5,  COMM_TAG_6,  COMM_TAG_7,  COMM_TAG_8
  use mpp_parameter_mod, only : COMM_TAG_9,  COMM_TAG_10, COMM_TAG_11, COMM_TAG_12
  use mpp_parameter_mod, only : COMM_TAG_13, COMM_TAG_14, COMM_TAG_15, COMM_TAG_16
  use mpp_parameter_mod, only : COMM_TAG_17, COMM_TAG_18, COMM_TAG_19, COMM_TAG_20
  use mpp_parameter_mod, only : MPP_FILL_INT,MPP_FILL_DOUBLE
  use mpp_data_mod,      only : stat, mpp_stack, ptr_stack, status, ptr_status, sync, ptr_sync
  use mpp_data_mod,      only : mpp_from_pe, ptr_from, remote_data_loc, ptr_remote
  use mpp_data_mod,      only : mpp_data_version=>version
  use platform_mod

implicit none
private

  !--- public parameters  -----------------------------------------------
  public :: MPP_VERBOSE, MPP_DEBUG, ALL_PES, ANY_PE, NULL_PE, NOTE, WARNING, FATAL
  public :: MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED, CLOCK_COMPONENT, CLOCK_SUBCOMPONENT
  public :: CLOCK_MODULE_DRIVER, CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA
  public :: MAXPES, EVENT_RECV, EVENT_SEND
  public :: COMM_TAG_1,  COMM_TAG_2,  COMM_TAG_3,  COMM_TAG_4
  public :: COMM_TAG_5,  COMM_TAG_6,  COMM_TAG_7,  COMM_TAG_8
  public :: COMM_TAG_9,  COMM_TAG_10, COMM_TAG_11, COMM_TAG_12
  public :: COMM_TAG_13, COMM_TAG_14, COMM_TAG_15, COMM_TAG_16
  public :: COMM_TAG_17, COMM_TAG_18, COMM_TAG_19, COMM_TAG_20
  public :: MPP_FILL_INT, MPP_FILL_DOUBLE, MPP_INFO_NULL, MPP_COMM_NULL
  public :: mpp_init_test_full_init, mpp_init_test_init_true_only, mpp_init_test_peset_allocated
  public :: mpp_init_test_clocks_init, mpp_init_test_datatype_list_init, mpp_init_test_logfile_init
  public :: mpp_init_test_read_namelist, mpp_init_test_etc_unit, mpp_init_test_requests_allocated

  !--- public interface from mpp_util.h ------------------------------
  public :: stdin, stdout, stderr, stdlog, warnlog, lowercase, uppercase, mpp_error, mpp_error_state
  public :: mpp_set_warn_level, mpp_sync, mpp_sync_self, mpp_pe
  public :: mpp_npes, mpp_root_pe, mpp_set_root_pe, mpp_declare_pelist
  public :: mpp_get_current_pelist, mpp_set_current_pelist, mpp_get_current_pelist_name
  public :: mpp_clock_id, mpp_clock_set_grain, mpp_record_timing_data, get_unit
  public :: read_ascii_file, read_input_nml, mpp_clock_begin, mpp_clock_end
  public :: get_ascii_file_num_lines, get_ascii_file_num_lines_and_length
  public :: mpp_record_time_start, mpp_record_time_end
  public :: mpp_commID, mpp_comm, inverse_permutation

  !--- public interface from mpp_comm.h ------------------------------
  public :: mpp_chksum, mpp_max, mpp_min, mpp_sum, mpp_transmit, mpp_send, mpp_recv
  public :: mpp_sum_ad
  public :: mpp_broadcast, mpp_init, mpp_exit
  public :: mpp_gather, mpp_scatter, mpp_alltoall
  public :: mpp_type, mpp_byte, mpp_type_create, mpp_type_free

  !*********************************************************************
  !
  !    public data type
  !
  !*********************************************************************
  !> Communication information for message passing libraries
  !!
  !> peset hold communicators as SHMEM-compatible triads (start, log2(stride), num)
  type :: communicator
     private
     character(len=32) :: name
     integer, pointer  :: list(:) =>NULL()
     integer           :: count
     integer           :: start, log2stride !< dummy variables when libMPI is defined.
     type(mpi_comm)    :: comm              !< MPI communicator for this PE set
     type(mpi_group)   :: group             !< MPI group for this PE set
  end type communicator

  !> Communication event profile
  type :: event
     private
     character(len=16)                         :: name
     integer(i8_kind), dimension(MAX_EVENTS)   :: ticks, bytes
     integer                                   :: calls
  end type event

  !> a clock contains an array of event profiles for a region
  type :: clock
     private
     character(len=32)    :: name
     integer(i8_kind)     :: hits
     integer(i8_kind)     :: tick
     integer(i8_kind)     :: total_ticks
     integer              :: peset_num
     logical              :: sync_on_begin, detailed
     integer              :: grain
     type(event), pointer :: events(:) =>NULL() !> if needed, allocate to MAX_EVENT_TYPES
     logical              :: is_on              !> initialize to false. set true when calling mpp_clock_begin
                                                !! set false when calling mpp_clock_end
  end type clock

  !> Summary of information from a clock run
  type :: Clock_Data_Summary
     private
     character(len=16)  :: name
     real(r8_kind)      :: msg_size_sums(MAX_BINS)
     real(r8_kind)      :: msg_time_sums(MAX_BINS)
     real(r8_kind)      :: total_data
     real(r8_kind)      :: total_time
     integer(i8_kind)   :: msg_size_cnts(MAX_BINS)
     integer(i8_kind)   :: total_cnts
  end type Clock_Data_Summary

  !> holds name and clock data for use in @ref mpp_util.h
  type :: Summary_Struct
     private
     character(len=16)         :: name
     type (Clock_Data_Summary) :: event(MAX_EVENT_TYPES)
  end type Summary_Struct

  !> Data types for generalized data transfer (e.g. MPI_Type)
  type :: mpp_type
     private
     integer :: counter !> Number of instances of this type
     integer :: ndims
     integer, allocatable :: sizes(:)
     integer, allocatable :: subsizes(:)
     integer, allocatable :: starts(:)
     type(mpi_datatype) :: etype   !> Elementary data type (e.g. MPI_BYTE)
     type(mpi_datatype) :: id      !> Identifier within message passing library (e.g. MPI)

     type(mpp_type), pointer :: prev => null()
     type(mpp_type), pointer :: next => null()
  end type mpp_type

  !> Persisent elements for linked list interaction
  type :: mpp_type_list
      private
      type(mpp_type), pointer :: head => null()
      type(mpp_type), pointer :: tail => null()
      integer :: length
  end type mpp_type_list

!***********************************************************************
!
!     public interface from mpp_util.h
!
!***********************************************************************
  !> @brief Error handler.
  !!
  !>    It is strongly recommended that all error exits pass through
  !!    <TT>mpp_error</TT> to assure the program fails cleanly. An individual
  !!    PE encountering a <TT>STOP</TT> statement, for instance, can cause the
  !!    program to hang. The use of the <TT>STOP</TT> statement is strongly
  !!    discouraged.
  !!
  !!    Calling mpp_error with no arguments produces an immediate error
  !!    exit, i.e:
  !!    <PRE>
  !!                    call mpp_error
  !!                    call mpp_error()
  !!    </PRE>
  !!    are equivalent.
  !!
  !!    The argument order
  !!    <PRE>
  !!                    call mpp_error( routine, errormsg, errortype )
  !!    </PRE>
  !!    is also provided to support legacy code. In this version of the
  !!    call, none of the arguments may be omitted.
  !!
  !!    The behaviour of <TT>mpp_error</TT> for a <TT>WARNING</TT> can be
  !!    controlled with an additional call <TT>mpp_set_warn_level</TT>.
  !!    <PRE>
  !!                    call mpp_set_warn_level(ERROR)
  !!    </PRE>
  !!    causes <TT>mpp_error</TT> to treat <TT>WARNING</TT>
  !!    exactly like <TT>FATAL</TT>.
  !!    <PRE>
  !!                    call mpp_set_warn_level(WARNING)
  !!    </PRE>
  !!    resets to the default behaviour described above.
  !!
  !!    <TT>mpp_error</TT> also has an internal error state which
  !!    maintains knowledge of whether a warning has been issued. This can be
  !!    used at startup in a subroutine that checks if the model has been
  !!    properly configured. You can generate a series of warnings using
  !!    <TT>mpp_error</TT>, and then check at the end if any warnings has been
  !!    issued using the function <TT>mpp_error_state()</TT>. If the value of
  !!    this is <TT>WARNING</TT>, at least one warning has been issued, and
  !!    the user can take appropriate action:
  !!
  !!    <PRE>
  !!                    if( ... )call mpp_error( WARNING, '...' )
  !!                    if( ... )call mpp_error( WARNING, '...' )
  !!                    if( ... )call mpp_error( WARNING, '...' )
  !!                    ...
  !!                    if( mpp_error_state().EQ.WARNING )call mpp_error( FATAL, '...' )
  !!    </PRE>
  !!  </DESCRIPTION>
  !! <br> Example usage:
  !! @code{.F90}
  !! call mpp_error( errortype, routine, errormsg )
  !! @endcode
  !! @param errortype
  !!    One of <TT>NOTE</TT>, <TT>WARNING</TT> or <TT>FATAL</TT>
  !!    (these definitions are acquired by use association).
  !!    <TT>NOTE</TT> writes <TT>errormsg</TT> to <TT>STDOUT</TT>.
  !!    <TT>WARNING</TT> writes <TT>errormsg</TT> to <TT>STDERR</TT>.
  !!    <TT>FATAL</TT> writes <TT>errormsg</TT> to <TT>STDERR</TT>,
  !!    and induces a clean error exit with a call stack traceback.
  !! @param routine Calling routine name
  !! @param errmsg Message to output
  !!  </IN>
  interface mpp_error
     module procedure mpp_error_basic
     module procedure mpp_error_mesg
     module procedure mpp_error_noargs
     module procedure mpp_error_is
     module procedure mpp_error_rs
     module procedure mpp_error_ia
     module procedure mpp_error_ra
     module procedure mpp_error_ia_ia
     module procedure mpp_error_ia_ra
     module procedure mpp_error_ra_ia
     module procedure mpp_error_ra_ra
     module procedure mpp_error_ia_is
     module procedure mpp_error_ia_rs
     module procedure mpp_error_ra_is
     module procedure mpp_error_ra_rs
     module procedure mpp_error_is_ia
     module procedure mpp_error_is_ra
     module procedure mpp_error_rs_ia
     module procedure mpp_error_rs_ra
     module procedure mpp_error_is_is
     module procedure mpp_error_is_rs
     module procedure mpp_error_rs_is
     module procedure mpp_error_rs_rs
  end interface

  !> Takes a given integer or real array and returns it as a string
  !> @param[in] array An array of integers or reals
  !> @returns string equivalent of given array

  interface array_to_char
     module procedure iarray_to_char
     module procedure rarray_to_char
  end interface

  !> Declare a pelist. The two flavors of this subroutine differ in the type
  !! of their comm/commID argument: mpp_declare_pelist_f08 expects a type(mpi_comm)
  !! as its comm argument, whereas mpp_declare_pelist_legacy expects an integer
  !! as its commID argument.
  interface mpp_declare_pelist
    module procedure mpp_declare_pelist_f08
    module procedure mpp_declare_pelist_legacy
  end interface

  !> Get the current pelist. The two flavors of this subroutine differ in the type
  !! of their comm/commID argument: mpp_get_current_pelist_f08 expects a type(mpi_comm)
  !! as its comm argument, whereas mpp_get_current_pelist_legacy expects an integer
  !! as its commID argument.
  interface mpp_get_current_pelist
    module procedure mpp_get_current_pelist_f08
    module procedure mpp_get_current_pelist_legacy
  end interface

!***********************************************************************
!
!    public interface from mpp_comm.h
!
!***********************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                             !
  !       ROUTINES TO INITIALIZE/FINALIZE MPP MODULE: mpp_init, mpp_exit        !
  !                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @fn mpp_mod::mpp_init::mpp_init( flags, localcomm, test_level)
!> @brief Initialize @ref mpp_mod
!!
!> Called to initialize the <TT>mpp_mod</TT> package. It is recommended
!! that this call be the first executed line in your program. It sets the
!! number of PEs assigned to this run (acquired from the command line, or
!! through the environment variable <TT>NPES</TT>), and associates an ID
!! number to each PE. These can be accessed by calling @ref mpp_npes and
!! @ref mpp_pe.
!! <br> Example usage:
!!
!!            call mpp_init( flags )
!!
!! @param flags
!!   <TT>flags</TT> can be set to <TT>MPP_VERBOSE</TT> to
!!   have <TT>mpp_mod</TT> keep you informed of what it's up to.
!! @param localcomm
!!   This is a type(mpi_comm) in mpp_init_f08, and an integer in mpp_init_legacy.
!!   This argument should only be used if MPI has previously been initialized by
!!   an external call to MPI_Init.
!! @param test_level
!!   Debugging flag to set amount of initialization tasks performed
  interface mpp_init
    module procedure mpp_init_f08
    module procedure mpp_init_legacy
  end interface

!> @fn mpp_mod::mpp_exit()
!> @brief Exit <TT>@ref mpp_mod</TT>.
!!
!> Called at the end of the run, or to re-initialize <TT>mpp_mod</TT>,
!! should you require that for some odd reason.
!!
!! This call implies synchronization across all PEs.
!!
!! <br>Example usage:
!!
!!            call mpp_exit()

  !#####################################################################

  !> @fn subroutine mpp_set_stack_size(n)
  !> @brief Allocate module internal workspace.
  !> @param Integer to set stack size to(in words)
  !> <TT>mpp_mod</TT> maintains a private internal array called
  !! <TT>mpp_stack</TT> for private workspace. This call sets the length,
  !! in words, of this array.
  !!
  !! The <TT>mpp_init</TT> call sets this
  !! workspace length to a default of 32768, and this call may be used if a
  !! longer workspace is needed.
  !!
  !! This call implies synchronization across all PEs.
  !!
  !! This workspace is symmetrically allocated, as required for
  !! efficient communication on SGI and Cray MPP systems. Since symmetric
  !! allocation must be performed by <I>all</I> PEs in a job, this call
  !! must also be called by all PEs, using the same value of
  !! <TT>n</TT>. Calling <TT>mpp_set_stack_size</TT> from a subset of PEs,
  !! or with unequal argument <TT>n</TT>, may cause the program to hang.
  !!
  !! If any MPP call using <TT>mpp_stack</TT> overflows the declared
  !! stack array, the program will abort with a message specifying the
  !! stack length that is required. Many users wonder why, if the required
  !! stack length can be computed, it cannot also be specified at that
  !! point. This cannot be automated because there is no way for the
  !! program to know if all PEs are present at that call, and with equal
  !! values of <TT>n</TT>. The program must be rerun by the user with the
  !! correct argument to <TT>mpp_set_stack_size</TT>, called at an
  !! appropriate point in the code where all PEs are known to be present.
  !!        @verbose call mpp_set_stack_size(n)
  !!
  public :: mpp_set_stack_size
  ! from mpp_util.h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!              DATA TRANSFER TYPES: mpp_type_create                           !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief Create a mpp_type variable
  !> @param[in] field A field of any numerical or logical type
  !> @param[in] array_of_subsizes Integer array of subsizes
  !> @param[in] array_of_starts Integer array of starts
  !> @param[out] dtype_out Output variable for created @ref mpp_type
  interface mpp_type_create
      module procedure mpp_type_create_int4
      module procedure mpp_type_create_int8
      module procedure mpp_type_create_real4
      module procedure mpp_type_create_real8
      module procedure mpp_type_create_cmplx4
      module procedure mpp_type_create_cmplx8
      module procedure mpp_type_create_logical4
      module procedure mpp_type_create_logical8
  end interface mpp_type_create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                             !
  !            GLOBAL REDUCTION ROUTINES: mpp_max, mpp_sum, mpp_min             !
  !                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief Reduction operations.
  !>    Find the max of scalar a from the PEs in pelist
  !!    result is also automatically broadcast to all PEs
  !!    @code{.F90}
  !!            call  mpp_max( a, pelist )
  !!    @endcode
  !> @param a <TT>real</TT> or <TT>integer</TT>, of 4-byte of 8-byte kind.
  !> @param pelist If <TT>pelist</TT> is omitted, the context is assumed to be the
  !!    current pelist. This call implies synchronization across the PEs in
  !!    <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
  interface mpp_max
     module procedure mpp_max_real8_0d
     module procedure mpp_max_real8_1d
     module procedure mpp_max_int8_0d
     module procedure mpp_max_int8_1d
     module procedure mpp_max_real4_0d
     module procedure mpp_max_real4_1d
     module procedure mpp_max_int4_0d
     module procedure mpp_max_int4_1d
  end interface

  !> @brief Reduction operations.
  !>    Find the min of scalar a from the PEs in pelist
  !!    result is also automatically broadcast to all PEs
  !!    @code{.F90}
  !!            call  mpp_min( a, pelist )
  !!    @endcode
  !> @param a <TT>real</TT> or <TT>integer</TT>, of 4-byte of 8-byte kind.
  !> @param pelist If <TT>pelist</TT> is omitted, the context is assumed to be the
  !!    current pelist. This call implies synchronization across the PEs in
  !!    <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
  interface mpp_min
     module procedure mpp_min_real8_0d
     module procedure mpp_min_real8_1d
     module procedure mpp_min_int8_0d
     module procedure mpp_min_int8_1d
     module procedure mpp_min_real4_0d
     module procedure mpp_min_real4_1d
     module procedure mpp_min_int4_0d
     module procedure mpp_min_int4_1d
  end interface


  !> @brief Reduction operation.
  !!
  !> <TT>MPP_TYPE_</TT> corresponds to any 4-byte and 8-byte variant of
  !! <TT>integer, real, complex</TT> variables, of rank 0 or 1. A
  !! contiguous block from a multi-dimensional array may be passed by its
  !! starting address and its length, as in <TT>f77</TT>.
  !!
  !! Library reduction operators are not required or guaranteed to be
  !! bit-reproducible. In any case, changing the processor count changes
  !! the data layout, and thus very likely the order of operations. For
  !! bit-reproducible sums of distributed arrays, consider using the
  !! <TT>mpp_global_sum</TT> routine provided by the
  !! @ref mpp_domains module.
  !!
  !! The <TT>bit_reproducible</TT> flag provided in earlier versions of
  !! this routine has been removed.
  !!
  !!
  !! If <TT>pelist</TT> is omitted, the context is assumed to be the
  !! current pelist. This call implies synchronization across the PEs in
  !! <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
  !! Example usage:
  !!            call mpp_sum( a, length, pelist )
  !!
  interface mpp_sum
     module procedure mpp_sum_int8
     module procedure mpp_sum_int8_scalar
     module procedure mpp_sum_int8_2d
     module procedure mpp_sum_int8_3d
     module procedure mpp_sum_int8_4d
     module procedure mpp_sum_int8_5d
     module procedure mpp_sum_real8
     module procedure mpp_sum_real8_scalar
     module procedure mpp_sum_real8_2d
     module procedure mpp_sum_real8_3d
     module procedure mpp_sum_real8_4d
     module procedure mpp_sum_real8_5d
# 634 "mpp/mpp.F90"
     module procedure mpp_sum_int4
     module procedure mpp_sum_int4_scalar
     module procedure mpp_sum_int4_2d
     module procedure mpp_sum_int4_3d
     module procedure mpp_sum_int4_4d
     module procedure mpp_sum_int4_5d
     module procedure mpp_sum_real4
     module procedure mpp_sum_real4_scalar
     module procedure mpp_sum_real4_2d
     module procedure mpp_sum_real4_3d
     module procedure mpp_sum_real4_4d
     module procedure mpp_sum_real4_5d
# 654 "mpp/mpp.F90"
  end interface

  !> Calculates sum of a given numerical array across pe's for adjoint domains
  interface mpp_sum_ad
     module procedure mpp_sum_int8_ad
     module procedure mpp_sum_int8_scalar_ad
     module procedure mpp_sum_int8_2d_ad
     module procedure mpp_sum_int8_3d_ad
     module procedure mpp_sum_int8_4d_ad
     module procedure mpp_sum_int8_5d_ad
     module procedure mpp_sum_real8_ad
     module procedure mpp_sum_real8_scalar_ad
     module procedure mpp_sum_real8_2d_ad
     module procedure mpp_sum_real8_3d_ad
     module procedure mpp_sum_real8_4d_ad
     module procedure mpp_sum_real8_5d_ad
# 678 "mpp/mpp.F90"
     module procedure mpp_sum_int4_ad
     module procedure mpp_sum_int4_scalar_ad
     module procedure mpp_sum_int4_2d_ad
     module procedure mpp_sum_int4_3d_ad
     module procedure mpp_sum_int4_4d_ad
     module procedure mpp_sum_int4_5d_ad
     module procedure mpp_sum_real4_ad
     module procedure mpp_sum_real4_scalar_ad
     module procedure mpp_sum_real4_2d_ad
     module procedure mpp_sum_real4_3d_ad
     module procedure mpp_sum_real4_4d_ad
     module procedure mpp_sum_real4_5d_ad
# 698 "mpp/mpp.F90"
  end interface

  !> @brief Gather data sent from pelist onto the root pe
  !! Wrapper for MPI_gather, can be used with and without indices
  !!
  !> @param sbuf MPP_TYPE_ data buffer to send
  !> @param rbuf MPP_TYPE_ data buffer to receive
  !> @param pelist integer(:) optional pelist to gather from, defaults to current
  !>
  !> <BR> Example usage:
  !!
  !!                    call mpp_gather(send_buffer,recv_buffer, pelist)
  !!                    call mpp_gather(is, ie, js, je, pelist, array_seg, data, is_root_pe)
  !!
  interface mpp_gather
     module procedure mpp_gather_logical4
     module procedure mpp_gatherv_logical4
     module procedure mpp_gather_logical_1d
     module procedure mpp_gather_int4
     module procedure mpp_gather_int8
     module procedure mpp_gatherv_int4
     module procedure mpp_gatherv_int8
     module procedure mpp_gather_int4_1d
     module procedure mpp_gather_int8_1d
     module procedure mpp_gather_real4
     module procedure mpp_gather_real8
     module procedure mpp_gatherv_real4
     module procedure mpp_gatherv_real8
     module procedure mpp_gather_real4_1d
     module procedure mpp_gather_real8_1d
     module procedure mpp_gather_logical_1dv
     module procedure mpp_gather_int4_1dv
     module procedure mpp_gather_int8_1dv
     module procedure mpp_gather_real4_1dv
     module procedure mpp_gather_real8_1dv
     module procedure mpp_gather_pelist_logical_2d
     module procedure mpp_gather_pelist_logical_gen_2d
     module procedure mpp_gather_pelist_logical_3d
     module procedure mpp_gather_pelist_logical_gen_3d
     module procedure mpp_gather_pelist_int4_2d
     module procedure mpp_gather_pelist_int4_gen_2d
     module procedure mpp_gather_pelist_int4_3d
     module procedure mpp_gather_pelist_int4_gen_3d
     module procedure mpp_gather_pelist_int8_2d
     module procedure mpp_gather_pelist_int8_gen_2d
     module procedure mpp_gather_pelist_int8_3d
     module procedure mpp_gather_pelist_int8_gen_3d
     module procedure mpp_gather_pelist_real4_2d
     module procedure mpp_gather_pelist_real4_gen_2d
     module procedure mpp_gather_pelist_real4_3d
     module procedure mpp_gather_pelist_real4_gen_3d
     module procedure mpp_gather_pelist_real8_2d
     module procedure mpp_gather_pelist_real8_gen_2d
     module procedure mpp_gather_pelist_real8_3d
     module procedure mpp_gather_pelist_real8_gen_3d
  end interface

  !> @brief Scatter (ie - is) * (je - js) contiguous elements of array data from the designated root pe
  !! into contigous members of array segment in each pe that is included in the pelist argument.
  !!
  !> @param is, ie integer start and end index of the first dimension of the segment array
  !> @param je, js integer start and end index of the second dimension of the segment array
  !> @param pelist integer(:) the PE list of target pes, needs to be monotonically increasing
  !> @param array_seg MPP_TYPE_ 2D array that the data is to be copied into
  !> @param data MPP_TYPE_ the source array
  !> @param is_root_pe logical true if calling from root pe
  !> @param ishift integer offsets specifying the first elelement in the data array
  !> @param nk integer size of third dimension for 3D calls
  !!
  !> <BR> Example usage:
  !!
  !!                    call mpp_scatter(is, ie, js, je, pelist, segment, data, .true.)
  !!
  interface mpp_scatter
     module procedure mpp_scatterv_int4
     module procedure mpp_scatter_pelist_int4_2d
     module procedure mpp_scatter_pelist_int4_gen_2d
     module procedure mpp_scatter_pelist_int4_3d
     module procedure mpp_scatter_pelist_int4_gen_3d
     module procedure mpp_scatterv_int8
     module procedure mpp_scatter_pelist_int8_2d
     module procedure mpp_scatter_pelist_int8_gen_2d
     module procedure mpp_scatter_pelist_int8_3d
     module procedure mpp_scatter_pelist_int8_gen_3d
     module procedure mpp_scatterv_real4
     module procedure mpp_scatter_pelist_real4_2d
     module procedure mpp_scatter_pelist_real4_gen_2d
     module procedure mpp_scatter_pelist_real4_3d
     module procedure mpp_scatter_pelist_real4_gen_3d
     module procedure mpp_scatterv_real8
     module procedure mpp_scatter_pelist_real8_2d
     module procedure mpp_scatter_pelist_real8_gen_2d
     module procedure mpp_scatter_pelist_real8_3d
     module procedure mpp_scatter_pelist_real8_gen_3d
  end interface

  !#####################################################################
  !> @brief Scatter a vector across all PEs
  !!
  !> Transpose the vector and PE index
  !! Wrapper for the MPI_alltoall function, includes more generic _V and _W
  !! versions if given displacements/data types
  !!
  !! Generic MPP_TYPE_ implentations:
  !! <li> @ref mpp_alltoall_ </li>
  !! <li> @ref mpp_alltoallv_ </li>
  !! <li> @ref mpp_alltoallw_ </li>
  !!
  interface mpp_alltoall
     module procedure mpp_alltoall_int4
     module procedure mpp_alltoall_int8
     module procedure mpp_alltoall_real4
     module procedure mpp_alltoall_real8






     module procedure mpp_alltoall_logical4
     module procedure mpp_alltoall_logical8
     module procedure mpp_alltoall_int4_v
     module procedure mpp_alltoall_int8_v
     module procedure mpp_alltoall_real4_v
     module procedure mpp_alltoall_real8_v






     module procedure mpp_alltoall_logical4_v
     module procedure mpp_alltoall_logical8_v
     module procedure mpp_alltoall_int4_w
     module procedure mpp_alltoall_int8_w
     module procedure mpp_alltoall_real4_w
     module procedure mpp_alltoall_real8_w






     module procedure mpp_alltoall_logical4_w
     module procedure mpp_alltoall_logical8_w
  end interface


  !#####################################################################
  !> @brief Basic message-passing call.
  !!
  !>    <TT>MPP_TYPE_</TT> corresponds to any 4-byte and 8-byte variant of
  !!    <TT>integer, real, complex, logical</TT> variables, of rank 0 or 1. A
  !!    contiguous block from a multi-dimensional array may be passed by its
  !!    starting address and its length, as in <TT>f77</TT>.
  !!
  !!    <TT>mpp_transmit</TT> is currently implemented as asynchronous
  !!    outward transmission and synchronous inward transmission. This follows
  !!    the behaviour of <TT>shmem_put</TT> and <TT>shmem_get</TT>. In MPI, it
  !!    is implemented as <TT>mpi_isend</TT> and <TT>mpi_recv</TT>. For most
  !!    applications, transmissions occur in pairs, and are here accomplished
  !!    in a single call.
  !!
  !!    The special PE designations <TT>NULL_PE</TT>,
  !!    <TT>ANY_PE</TT> and <TT>ALL_PES</TT> are provided by use
  !!    association.
  !!
  !!    <TT>NULL_PE</TT>: is used to disable one of the pair of
  !!    transmissions.<BR/>
  !!    <TT>ANY_PE</TT>: is used for unspecific remote
  !!    destination. (Please note that <TT>put_pe=ANY_PE</TT> has no meaning
  !!    in the MPI context, though it is available in the SHMEM invocation. If
  !!    portability is a concern, it is best avoided).<BR/>
  !!    <TT>ALL_PES</TT>: is used for broadcast operations.
  !!
  !!    It is recommended that
  !!    @ref mpp_broadcast be used for
  !!    broadcasts.
  !!
  !!    The following example illustrates the use of
  !!    <TT>NULL_PE</TT> and <TT>ALL_PES</TT>:
  !!
  !!    <PRE>
  !!    real, dimension(n) :: a
  !!    if( pe.EQ.0 )then
  !!        do p = 1,npes-1
  !!           call mpp_transmit( a, n, p, a, n, NULL_PE )
  !!        end do
  !!    else
  !!        call mpp_transmit( a, n, NULL_PE, a, n, 0 )
  !!    end if
  !!
  !!    call mpp_transmit( a, n, ALL_PES, a, n, 0 )
  !!    </PRE>
  !!
  !!    The do loop and the broadcast operation above are equivalent.
  !!
  !!    Two overloaded calls <TT>mpp_send</TT> and
  !!     <TT>mpp_recv</TT> have also been
  !!    provided. <TT>mpp_send</TT> calls <TT>mpp_transmit</TT>
  !!    with <TT>get_pe=NULL_PE</TT>. <TT>mpp_recv</TT> calls
  !!    <TT>mpp_transmit</TT> with <TT>put_pe=NULL_PE</TT>. Thus
  !!    the do loop above could be written more succinctly:
  !!
  !!    <PRE>
  !!    if( pe.EQ.0 )then
  !!        do p = 1,npes-1
  !!           call mpp_send( a, n, p )
  !!        end do
  !!    else
  !!        call mpp_recv( a, n, 0 )
  !!    end if
  !!    </PRE>
  !! <br>Example call:
  !! @code{.F90}
  !!    call mpp_transmit( put_data, put_len, put_pe, get_data, get_len, get_pe )
  !! @endcode
  interface mpp_transmit
     module procedure mpp_transmit_real8
     module procedure mpp_transmit_real8_scalar
     module procedure mpp_transmit_real8_2d
     module procedure mpp_transmit_real8_3d
     module procedure mpp_transmit_real8_4d
     module procedure mpp_transmit_real8_5d
# 930 "mpp/mpp.F90"
     module procedure mpp_transmit_int8
     module procedure mpp_transmit_int8_scalar
     module procedure mpp_transmit_int8_2d
     module procedure mpp_transmit_int8_3d
     module procedure mpp_transmit_int8_4d
     module procedure mpp_transmit_int8_5d
     module procedure mpp_transmit_logical8
     module procedure mpp_transmit_logical8_scalar
     module procedure mpp_transmit_logical8_2d
     module procedure mpp_transmit_logical8_3d
     module procedure mpp_transmit_logical8_4d
     module procedure mpp_transmit_logical8_5d

     module procedure mpp_transmit_real4
     module procedure mpp_transmit_real4_scalar
     module procedure mpp_transmit_real4_2d
     module procedure mpp_transmit_real4_3d
     module procedure mpp_transmit_real4_4d
     module procedure mpp_transmit_real4_5d

# 958 "mpp/mpp.F90"
     module procedure mpp_transmit_int4
     module procedure mpp_transmit_int4_scalar
     module procedure mpp_transmit_int4_2d
     module procedure mpp_transmit_int4_3d
     module procedure mpp_transmit_int4_4d
     module procedure mpp_transmit_int4_5d
     module procedure mpp_transmit_logical4
     module procedure mpp_transmit_logical4_scalar
     module procedure mpp_transmit_logical4_2d
     module procedure mpp_transmit_logical4_3d
     module procedure mpp_transmit_logical4_4d
     module procedure mpp_transmit_logical4_5d
  end interface
  !> @brief Receive data from another PE
  !!
  !> @param[out] get_data scalar or array to get written with received data
  !> @param get_len size of array to recv from get_data
  !> @param from_pe PE number to receive from
  !> @param block true for blocking, false for non-blocking. Defaults to true
  !> @param tag communication tag
  !> @param[out] request MPI request handle
  interface mpp_recv
     module procedure mpp_recv_real8
     module procedure mpp_recv_real8_scalar
     module procedure mpp_recv_real8_2d
     module procedure mpp_recv_real8_3d
     module procedure mpp_recv_real8_4d
     module procedure mpp_recv_real8_5d
# 994 "mpp/mpp.F90"
     module procedure mpp_recv_int8
     module procedure mpp_recv_int8_scalar
     module procedure mpp_recv_int8_2d
     module procedure mpp_recv_int8_3d
     module procedure mpp_recv_int8_4d
     module procedure mpp_recv_int8_5d
     module procedure mpp_recv_logical8
     module procedure mpp_recv_logical8_scalar
     module procedure mpp_recv_logical8_2d
     module procedure mpp_recv_logical8_3d
     module procedure mpp_recv_logical8_4d
     module procedure mpp_recv_logical8_5d

     module procedure mpp_recv_real4
     module procedure mpp_recv_real4_scalar
     module procedure mpp_recv_real4_2d
     module procedure mpp_recv_real4_3d
     module procedure mpp_recv_real4_4d
     module procedure mpp_recv_real4_5d

# 1022 "mpp/mpp.F90"
     module procedure mpp_recv_int4
     module procedure mpp_recv_int4_scalar
     module procedure mpp_recv_int4_2d
     module procedure mpp_recv_int4_3d
     module procedure mpp_recv_int4_4d
     module procedure mpp_recv_int4_5d
     module procedure mpp_recv_logical4
     module procedure mpp_recv_logical4_scalar
     module procedure mpp_recv_logical4_2d
     module procedure mpp_recv_logical4_3d
     module procedure mpp_recv_logical4_4d
     module procedure mpp_recv_logical4_5d
  end interface
  !> Send data to a receiving PE.
  !!
  !> @param put_data scalar or array to get sent to a receiving PE
  !> @param put_len size of data to send from put_data
  !> @param to_pe PE number to send to
  !> @param block true for blocking, false for non-blocking. Defaults to true
  !> @param tag communication tag
  !> @param[out] request MPI request handle
  !! <br> Example usage:
  !! @code{.F90} call mpp_send(data, ie, pe) @endcode
  interface mpp_send
     module procedure mpp_send_real8
     module procedure mpp_send_real8_scalar
     module procedure mpp_send_real8_2d
     module procedure mpp_send_real8_3d
     module procedure mpp_send_real8_4d
     module procedure mpp_send_real8_5d
# 1060 "mpp/mpp.F90"
     module procedure mpp_send_int8
     module procedure mpp_send_int8_scalar
     module procedure mpp_send_int8_2d
     module procedure mpp_send_int8_3d
     module procedure mpp_send_int8_4d
     module procedure mpp_send_int8_5d
     module procedure mpp_send_logical8
     module procedure mpp_send_logical8_scalar
     module procedure mpp_send_logical8_2d
     module procedure mpp_send_logical8_3d
     module procedure mpp_send_logical8_4d
     module procedure mpp_send_logical8_5d

     module procedure mpp_send_real4
     module procedure mpp_send_real4_scalar
     module procedure mpp_send_real4_2d
     module procedure mpp_send_real4_3d
     module procedure mpp_send_real4_4d
     module procedure mpp_send_real4_5d

# 1088 "mpp/mpp.F90"
     module procedure mpp_send_int4
     module procedure mpp_send_int4_scalar
     module procedure mpp_send_int4_2d
     module procedure mpp_send_int4_3d
     module procedure mpp_send_int4_4d
     module procedure mpp_send_int4_5d
     module procedure mpp_send_logical4
     module procedure mpp_send_logical4_scalar
     module procedure mpp_send_logical4_2d
     module procedure mpp_send_logical4_3d
     module procedure mpp_send_logical4_4d
     module procedure mpp_send_logical4_5d
  end interface


  !> @brief Perform parallel broadcasts
  !!
  !> The <TT>mpp_broadcast</TT> call has been added because the original
  !! syntax (using <TT>ALL_PES</TT> in <TT>mpp_transmit</TT>) did not
  !! support a broadcast across a pelist.
  !!
  !! <TT>MPP_TYPE_</TT> corresponds to any 4-byte and 8-byte variant of
  !! <TT>integer, real, complex, logical</TT> variables, of rank 0 or 1. A
  !! contiguous block from a multi-dimensional array may be passed by its
  !! starting address and its length, as in <TT>f77</TT>.
  !!
  !! Global broadcasts through the <TT>ALL_PES</TT> argument to
  !! @ref mpp_transmit are still provided for
  !! backward-compatibility.
  !!
  !! If <TT>pelist</TT> is omitted, the context is assumed to be the
  !! current pelist. <TT>from_pe</TT> must belong to the current
  !! pelist. This call implies synchronization across the PEs in
  !! <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
  !!
  !! <br>Example usage:
  !!
  !!            call mpp_broadcast( data, length, from_pe, pelist )
  !!
  !> @param[inout] data Data to broadcast
  !> @param length Length of data to broadcast
  !> @param from_pe PE to send the data from
  !> @param pelist List of PE's to broadcast across, if not provided uses current list
  interface mpp_broadcast
     module procedure mpp_broadcast_char
     module procedure mpp_broadcast_real8
     module procedure mpp_broadcast_real8_scalar
     module procedure mpp_broadcast_real8_2d
     module procedure mpp_broadcast_real8_3d
     module procedure mpp_broadcast_real8_4d
     module procedure mpp_broadcast_real8_5d
# 1147 "mpp/mpp.F90"
     module procedure mpp_broadcast_int8
     module procedure mpp_broadcast_int8_scalar
     module procedure mpp_broadcast_int8_2d
     module procedure mpp_broadcast_int8_3d
     module procedure mpp_broadcast_int8_4d
     module procedure mpp_broadcast_int8_5d
     module procedure mpp_broadcast_logical8
     module procedure mpp_broadcast_logical8_scalar
     module procedure mpp_broadcast_logical8_2d
     module procedure mpp_broadcast_logical8_3d
     module procedure mpp_broadcast_logical8_4d
     module procedure mpp_broadcast_logical8_5d

     module procedure mpp_broadcast_real4
     module procedure mpp_broadcast_real4_scalar
     module procedure mpp_broadcast_real4_2d
     module procedure mpp_broadcast_real4_3d
     module procedure mpp_broadcast_real4_4d
     module procedure mpp_broadcast_real4_5d

# 1175 "mpp/mpp.F90"
     module procedure mpp_broadcast_int4
     module procedure mpp_broadcast_int4_scalar
     module procedure mpp_broadcast_int4_2d
     module procedure mpp_broadcast_int4_3d
     module procedure mpp_broadcast_int4_4d
     module procedure mpp_broadcast_int4_5d
     module procedure mpp_broadcast_logical4
     module procedure mpp_broadcast_logical4_scalar
     module procedure mpp_broadcast_logical4_2d
     module procedure mpp_broadcast_logical4_3d
     module procedure mpp_broadcast_logical4_4d
     module procedure mpp_broadcast_logical4_5d
  end interface

  !#####################################################################

  !> @brief Calculate parallel checksums
  !!
  !> \e mpp_chksum is a parallel checksum routine that returns an
  !! identical answer for the same array irrespective of how it has been
  !! partitioned across processors. \e int_kind is the KIND
  !! parameter corresponding to long integers (see discussion on
  !! OS-dependent preprocessor directives) defined in
  !! the file platform.F90. \e MPP_TYPE_ corresponds to any
  !! 4-byte and 8-byte variant of \e integer, \e real, \e complex, \e logical
  !! variables, of rank 0 to 5.
  !!
  !! Integer checksums on FP data use the F90 <TT>TRANSFER()</TT>
  !! intrinsic.
  !!
  !! This provides identical results on a single-processor job, and to perform
  !! serial checksums on a single processor of a parallel job, you only
  !! need to use the optional <TT>pelist</TT> argument.
  !! <PRE>
  !! use mpp_mod
  !! integer :: pe, chksum
  !! real :: a(:)
  !! pe = mpp_pe()
  !! chksum = mpp_chksum( a, (/pe/) )
  !! </PRE>
  !!
  !! The additional functionality of <TT>mpp_chksum</TT> over
  !! serial checksums is to compute the checksum across the PEs in
  !! <TT>pelist</TT>. The answer is guaranteed to be the same for
  !! the same distributed array irrespective of how it has been
  !! partitioned.
  !!
  !! If <TT>pelist</TT> is omitted, the context is assumed to be the
  !! current pelist. This call implies synchronization across the PEs in
  !! <TT>pelist</TT>, or the current pelist if <TT>pelist</TT> is absent.
  !! <br> Example usage:
  !!
  !!            mpp_chksum( var, pelist )
  !!
  !! @param var Data to calculate checksum of
  !! @param pelist Optional list of PE's to include in checksum calculation if not using
  !! current pelist
  !! @return Parallel checksum of var across given or implicit pelist
  !!
  !! Generic MPP_TYPE_ implentations:
  !! <li> @ref mpp_chksum_</li>
  !! <li> @ref mpp_chksum_int_</li>
  !! <li> @ref mpp_chksum_int_rmask_</li>
  !!
  interface mpp_chksum
     module procedure mpp_chksum_i8_1d
     module procedure mpp_chksum_i8_2d
     module procedure mpp_chksum_i8_3d
     module procedure mpp_chksum_i8_4d
     module procedure mpp_chksum_i8_5d
     module procedure mpp_chksum_i8_1d_rmask
     module procedure mpp_chksum_i8_2d_rmask
     module procedure mpp_chksum_i8_3d_rmask
     module procedure mpp_chksum_i8_4d_rmask
     module procedure mpp_chksum_i8_5d_rmask

     module procedure mpp_chksum_i4_1d
     module procedure mpp_chksum_i4_2d
     module procedure mpp_chksum_i4_3d
     module procedure mpp_chksum_i4_4d
     module procedure mpp_chksum_i4_5d
     module procedure mpp_chksum_i4_1d_rmask
     module procedure mpp_chksum_i4_2d_rmask
     module procedure mpp_chksum_i4_3d_rmask
     module procedure mpp_chksum_i4_4d_rmask
     module procedure mpp_chksum_i4_5d_rmask

     module procedure mpp_chksum_r8_0d
     module procedure mpp_chksum_r8_1d
     module procedure mpp_chksum_r8_2d
     module procedure mpp_chksum_r8_3d
     module procedure mpp_chksum_r8_4d
     module procedure mpp_chksum_r8_5d

     module procedure mpp_chksum_r4_0d
     module procedure mpp_chksum_r4_1d
     module procedure mpp_chksum_r4_2d
     module procedure mpp_chksum_r4_3d
     module procedure mpp_chksum_r4_4d
     module procedure mpp_chksum_r4_5d
# 1291 "mpp/mpp.F90"
  end interface

!***********************************************************************
!
!            module variables
!
!***********************************************************************
  integer, parameter   :: PESET_MAX = 10000
  integer              :: current_peset_max = 32
  type(communicator), allocatable :: peset(:) !< Will be allocated starting from 0, 0 is a dummy used
                                              !! to hold single-PE "self" communicator
  logical              :: module_is_initialized = .false.
  logical              :: debug = .false.
  integer              :: npes=1, root_pe=0, pe=0
  integer(i8_kind)     :: tick, ticks_per_sec, max_ticks, start_tick, end_tick, tick0=0
  type(mpi_comm)       :: mpp_comm_private
  logical              :: first_call_system_clock_mpi=.TRUE.
  real(r8_kind)        :: mpi_count0=0  !< use to prevent integer overflow
  real(r8_kind)        :: mpi_tick_rate=0.d0  !< clock rate for mpi_wtick()
  logical              :: mpp_record_timing_data=.TRUE.
  type(clock),save     :: clocks(MAX_CLOCKS)
  integer              :: log_unit, etc_unit
  integer              :: warn_unit !< unit number of the warning log
  character(len=32), parameter    :: configfile='logfile'
  character(len=32), parameter    :: warnfile='warnfile' !< base name for warninglog (appends ".<PE>.out")
  integer              :: peset_num=0, current_peset_num=0
  integer              :: world_peset_num                  !<the world communicator
  integer              :: error
  integer              :: clock_num=0, num_clock_ids=0,current_clock=0, previous_clock(MAX_CLOCKS)=0
  real                 :: tick_rate

  type(mpp_type_list)    :: datatypes
  type(mpp_type), target :: mpp_byte

  integer              :: cur_send_request = 0
  integer              :: cur_recv_request = 0
  type(mpi_request), allocatable :: request_send(:)
  type(mpi_request), allocatable :: request_recv(:)
  type(mpi_datatype), allocatable :: type_recv(:)
  integer, allocatable :: size_recv(:)
! if you want to save the non-root PE information uncomment out the following line
! and comment out the assigment of etcfile to '/dev/null'



  character(len=32)    :: etcfile='/dev/null'


!> Use the intrinsics in iso_fortran_env
  integer :: in_unit=INPUT_UNIT, out_unit=OUTPUT_UNIT, err_unit=ERROR_UNIT
  integer :: stdout_unit

  !--- variables used in mpp_util.h
  type(Summary_Struct) :: clock_summary(MAX_CLOCKS)
  logical              :: warnings_are_fatal = .FALSE.
  integer              :: error_state=0
  integer              :: clock_grain=CLOCK_LOOP-1

  !--- variables used in mpp_comm.h
  integer            :: clock0    !<measures total runtime from mpp_init to mpp_exit
  integer            :: mpp_stack_size=0, mpp_stack_hwm=0
  logical            :: verbose=.FALSE.

  integer :: get_len_nocomm = 0 !< needed for mpp_transmit_nocomm.h

  !--- variables used in mpp_comm_mpi.inc
  integer, parameter :: mpp_init_test_full_init = -1
  integer, parameter :: mpp_init_test_init_true_only = 0
  integer, parameter :: mpp_init_test_peset_allocated = 1
  integer, parameter :: mpp_init_test_clocks_init = 2
  integer, parameter :: mpp_init_test_datatype_list_init = 3
  integer, parameter :: mpp_init_test_logfile_init = 4
  integer, parameter :: mpp_init_test_read_namelist = 5
  integer, parameter :: mpp_init_test_etc_unit = 6
  integer, parameter :: mpp_init_test_requests_allocated = 7

  !> MPP_INFO_NULL acts as an analagous mpp-macro for MPI_INFO_NULL to share with fms2_io NetCDF4
  !! mpi-io. Intel MPI and MPICH provide a value of 469762048. OpenMPI provides a value of 0.
  integer, parameter ::  MPP_INFO_NULL = MPI_INFO_NULL%mpi_val

  !> MPP_COMM_NULL acts as an analagous mpp-macro for MPI_COMM_NULL to share with fms2_io NetCDF4
  !! mpi-io. Intel MPI and MPICH provide a value of 67108864. OpenMPI provides a value of 0.
  integer, parameter ::  MPP_COMM_NULL = MPI_COMM_NULL%mpi_val

!***********************************************************************
!  variables needed for subroutine read_input_nml (include/mpp_util.inc)
!
! public variable needed for reading input nml file from an internal file
  character(len=:), dimension(:), allocatable, target, public :: input_nml_file
  logical :: read_ascii_file_on = .FALSE.
!***********************************************************************

! Include variable "version" to be written to log file.

# 1 "./include/file_version.h" 1
! -*-f90-*-
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************




  character(len=*), parameter :: version = 'unknown'
# 1385 "mpp/mpp.F90" 2
  public version

  integer, parameter :: MAX_REQUEST_MIN  = 10000
  integer            :: request_multiply = 20

  logical :: etc_unit_is_stderr = .false.
  integer :: max_request = 0
  logical :: sync_all_clocks = .false.
  namelist /mpp_nml/ etc_unit_is_stderr, request_multiply, mpp_record_timing_data, sync_all_clocks

  contains

# 1 "mpp/include/system_clock.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

# 49 "mpp/include/system_clock.fh"
subroutine system_clock_default( count, count_rate, count_max )
!mimics F90 system_clock_default intrinsic
      integer(i8_kind), optional :: count, count_rate, count_max
!count must return a number between 0 and count_max
      integer                      :: count_int, count_rate_int, count_max_int
      call system_clock( count_int, count_rate_int, count_max_int)
      if( PRESENT(count) )      count      = count_int
      if( PRESENT(count_rate) ) count_rate = count_rate_int
      if( PRESENT(count_max) )  count_max  = count_max_int
      return
    end subroutine system_clock_default
# 1397 "mpp/mpp.F90" 2

# 1 "mpp/include/mpp_util.inc" 1
! -*-f90-*-


!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************





# 1 "mpp/include/mpp_util_nocomm.inc" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> @file
!> @brief Utility routines for parallelization, non-mpi version

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!         MISCELLANEOUS UTILITIES: mpp_error                                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mpp_error_basic( errortype, errormsg )
  !a very basic error handler
  !uses ABORT and FLUSH calls, may need to use cpp to rename
  integer,                    intent(in) :: errortype
  character(len=*), intent(in), optional :: errormsg
  character(len=512)                     :: text
  logical                                :: opened
  integer                                :: istat, errunit, outunit

  if( .NOT.module_is_initialized )call ABORT()

  select case( errortype )
  case(NOTE)
     text = 'NOTE'         !just FYI
  case(WARNING)
     text = 'WARNING'      !probable error
  case(FATAL)
     text = 'FATAL'        !fatal error
  case default
     text = 'WARNING: non-existent errortype (must be NOTE|WARNING|FATAL)'
  end select

  if( npes.GT.1 )write( text,'(a,i5)' )trim(text)//' from PE', pe   !this is the mpp part
  if( PRESENT(errormsg) )text = trim(text)//': '//trim(errormsg)

  errunit = stderr()
  outunit = stdout()

  select case( errortype )
  case(NOTE)
     write( outunit,'(a)' )trim(text)
  case default
     write( errunit,'(/a/)' )trim(text)
     write( outunit,'(/a/)' )trim(text)
     if( errortype.EQ.FATAL .OR. warnings_are_fatal )then
        FLUSH(outunit)
        call ABORT() !automatically calls traceback on Cray systems
     end if
  end select

  error_state = errortype
  return
end subroutine mpp_error_basic

!#####################################################################
!> Makes a PE set out of a PE list. A PE list is an ordered list of PEs
!! a PE set is a triad (start,log2stride,size) for SHMEM, an a communicator for MPI
!! if stride is non-uniform or not a power of 2,
!! will return error (not required for MPI but enforced for uniformity)
function get_peset(pelist)
  integer                       :: get_peset
  integer, intent(in), optional :: pelist(:)

  if( .NOT.PRESENT(pelist) )then !set it to current_peset_num
     get_peset = current_peset_num; return
  end if

  get_peset = 0

  return

end function get_peset

!#######################################################################
!> synchronize PEs in list
subroutine mpp_sync( pelist, check )
  integer, intent(in), optional :: pelist(:)
  integer, intent(in), optional :: check

  return
end subroutine mpp_sync

!#######################################################################
!> This is to check if current PE's outstanding puts are complete
!! but we can't use shmem_fence because we are actually waiting for
!! a remote PE to complete its get
subroutine mpp_sync_self( pelist, check, request, msg_size, msg_type )
  integer, intent(in), optional :: pelist(:)
  integer, intent(in), optional :: check
  type(mpi_request), intent(inout), optional :: request(:)
  integer, intent(in), optional :: msg_size(:)
  type(mpi_datatype), intent(in), optional :: msg_type(:)


  return
end subroutine mpp_sync_self
# 26 "mpp/include/mpp_util.inc" 2


  !> @brief This function returns the current standard fortran unit numbers for input.
  function stdin()
    integer :: stdin
    stdin = in_unit
    return
  end function stdin

  !> @brief This function returns the current  standard fortran unit numbers for output.
  function stdout()
    integer :: stdout
    stdout = out_unit
    if( pe.NE.root_pe )stdout = stdlog()
    return
  end function stdout

  !> @brief This function returns the current standard fortran unit numbers for error messages.
  function stderr()
    integer :: stderr
    stderr = err_unit
    return
  end function stderr

  !> @brief This function returns the current  standard fortran unit numbers for log messages.
  !!    Log messages, by convention, are written to the file <TT>logfile.out</TT>.
  function stdlog()
    integer :: stdlog
    logical :: opened
    character(len=11) :: this_pe
!$  logical           :: omp_in_parallel
!$  integer           :: omp_get_num_threads
!$  integer           :: errunit


!NOTES: We can not use mpp_error to handle the error because mpp_error
!       will call stdout and stdout will call stdlog for non-root-pe.
!       This will be a cicular call.

!$  if( omp_in_parallel() .and. (omp_get_num_threads() > 1) ) then
!$OMP single
!$      errunit = stderr()
!$      write( errunit,'(/a/)' ) 'FATAL: STDLOG: is called inside a OMP parallel region'



!$      call ABORT()

!$OMP end single
!$  endif

    if( pe.EQ.root_pe )then
       write(this_pe,'(a,i6.6,a)') '.',pe,'.out'
       inquire( file=trim(configfile)//this_pe, opened=opened )
       if( opened )then
          FLUSH(log_unit)
       else
          open(newunit=log_unit, status='UNKNOWN', file=trim(configfile)//this_pe, position='APPEND', err=10 )
       end if
       stdlog = log_unit
    else
       inquire(unit=etc_unit, opened=opened )
       if( opened )then
          FLUSH(etc_unit)
       else
          open(newunit=etc_unit, status='UNKNOWN', file=trim(etcfile), position='APPEND', err=11 )
       end if
       stdlog = etc_unit
    end if
    return
10  call mpp_error( FATAL, 'STDLOG: unable to open '//trim(configfile)//this_pe//'.' )
11  call mpp_error( FATAL, 'STDLOG: unable to open '//trim(etcfile)//'.' )
  end function stdlog

  !#####################################################################
  subroutine mpp_init_logfile()
  integer :: p
  logical :: exist
  character(len=11) :: this_pe
  if( pe.EQ.root_pe )then
     do p=0,npes-1
       write(this_pe,'(a,i6.6,a)') '.',p,'.out'
       inquire( file=trim(configfile)//this_pe, exist=exist )
       if(exist)then
         open(newunit=log_unit, file=trim(configfile)//this_pe, status='REPLACE' )
         close(log_unit)
       endif
     end do
  end if
  end subroutine mpp_init_logfile

  !> Opens the warning log file, called during mpp_init
  subroutine mpp_init_warninglog()
  logical :: exist
  character(len=11) :: this_pe
  if( pe.EQ.root_pe )then
    write(this_pe,'(a,i6.6,a)') '.',pe,'.out'
    inquire( file=trim(warnfile)//this_pe, exist=exist )
    if(exist)then
      open(newunit=warn_unit, file=trim(warnfile)//this_pe, status='REPLACE' )
    else
      open(newunit=warn_unit, file=trim(warnfile)//this_pe, status='NEW' )
    endif
  end if
  end subroutine mpp_init_warninglog

  !> @brief This function returns unit number for the warning log
  !! if on the root pe, otherwise returns the etc_unit value (usually /dev/null)
  function warnlog()
    integer :: warnlog
    if(.not. module_is_initialized) call mpp_error(FATAL, "mpp_mod: warnlog cannot be called before mpp_init")
    if(root_pe .eq. pe) then
      warnlog = warn_unit
    else
      warnlog = etc_unit
    endif
    return
  end function warnlog

  !#####################################################################
  subroutine mpp_set_warn_level(flag)
    integer, intent(in) :: flag

    if( flag.EQ.WARNING )then
       warnings_are_fatal = .FALSE.
    else if( flag.EQ.FATAL )then
       warnings_are_fatal = .TRUE.
    else
       call mpp_error( FATAL, 'MPP_SET_WARN_LEVEL: warning flag must be set to WARNING or FATAL.' )
    end if
    return
  end subroutine mpp_set_warn_level

  !#####################################################################
  function mpp_error_state()
    integer :: mpp_error_state
    mpp_error_state = error_state
    return
  end function mpp_error_state

!#####################################################################
!> @brief overloads to mpp_error_basic, support for error_mesg routine in FMS
subroutine mpp_error_mesg( routine, errormsg, errortype )
  character(len=*), intent(in) :: routine, errormsg
  integer,          intent(in) :: errortype

  call mpp_error( errortype, trim(routine)//': '//trim(errormsg) )
  return
end subroutine mpp_error_mesg

!#####################################################################
subroutine mpp_error_noargs()
  call mpp_error(FATAL)
end subroutine mpp_error_noargs

!#####################################################################
subroutine mpp_error_Is(errortype, errormsg1, mpp_ival, errormsg2)
  integer,          intent(in) :: errortype
  INTEGER,          intent(in) :: mpp_ival
  character(len=*), intent(in) :: errormsg1
  character(len=*),      intent(in), optional :: errormsg2
  call mpp_error( errortype, errormsg1, (/mpp_ival/), errormsg2)
end subroutine mpp_error_Is
!#####################################################################
subroutine mpp_error_Rs(errortype, errormsg1, mpp_rval, errormsg2)
  integer,          intent(in) :: errortype
  REAL,             intent(in) :: mpp_rval
  character(len=*), intent(in) :: errormsg1
  character(len=*),      intent(in), optional :: errormsg2
  call mpp_error( errortype, errormsg1, (/mpp_rval/), errormsg2)
end subroutine mpp_error_Rs
!#####################################################################
subroutine mpp_error_Ia(errortype, errormsg1, array, errormsg2)
  integer,               intent(in) :: errortype
  INTEGER, dimension(:), intent(in) :: array
  character(len=*),      intent(in) :: errormsg1
  character(len=*),      intent(in), optional :: errormsg2
  character(len=512) :: string

  string = errormsg1//trim(array_to_char(array))
  if(present(errormsg2)) string = trim(string)//errormsg2
  call mpp_error_basic( errortype, trim(string))

end subroutine mpp_error_Ia

!#####################################################################
subroutine mpp_error_Ra(errortype, errormsg1, array, errormsg2)
  integer,            intent(in) :: errortype
  REAL, dimension(:), intent(in) :: array
  character(len=*),      intent(in) :: errormsg1
  character(len=*),   intent(in), optional :: errormsg2
  character(len=512) :: string

  string = errormsg1//trim(array_to_char(array))
  if(present(errormsg2)) string = trim(string)//errormsg2
  call mpp_error_basic( errortype, trim(string))

end subroutine mpp_error_Ra

!#####################################################################




# 1 "mpp/include/mpp_error_a_a.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_error_ia_ia(errortype, errormsg1, array1, errormsg2, array2, errormsg3)
  integer,            intent(in) :: errortype
  integer, dimension(:), intent(in) :: array1
  integer, dimension(:), intent(in) :: array2
  character(len=*),      intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3
  character(len=512) :: string

  string = errormsg1//trim(array_to_char(array1))
  string = trim(string)//errormsg2//trim(array_to_char(array2))
  if(present(errormsg3)) string = trim(string)//errormsg3
  call mpp_error_basic( errortype, trim(string))

end subroutine mpp_error_ia_ia
# 230 "mpp/include/mpp_util.inc" 2



!#####################################################################




# 1 "mpp/include/mpp_error_a_a.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_error_ia_ra(errortype, errormsg1, array1, errormsg2, array2, errormsg3)
  integer,            intent(in) :: errortype
  integer, dimension(:), intent(in) :: array1
  real, dimension(:), intent(in) :: array2
  character(len=*),      intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3
  character(len=512) :: string

  string = errormsg1//trim(array_to_char(array1))
  string = trim(string)//errormsg2//trim(array_to_char(array2))
  if(present(errormsg3)) string = trim(string)//errormsg3
  call mpp_error_basic( errortype, trim(string))

end subroutine mpp_error_ia_ra
# 238 "mpp/include/mpp_util.inc" 2



!#####################################################################




# 1 "mpp/include/mpp_error_a_a.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_error_ra_ia(errortype, errormsg1, array1, errormsg2, array2, errormsg3)
  integer,            intent(in) :: errortype
  real, dimension(:), intent(in) :: array1
  integer, dimension(:), intent(in) :: array2
  character(len=*),      intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3
  character(len=512) :: string

  string = errormsg1//trim(array_to_char(array1))
  string = trim(string)//errormsg2//trim(array_to_char(array2))
  if(present(errormsg3)) string = trim(string)//errormsg3
  call mpp_error_basic( errortype, trim(string))

end subroutine mpp_error_ra_ia
# 246 "mpp/include/mpp_util.inc" 2



!#####################################################################




# 1 "mpp/include/mpp_error_a_a.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_error_ra_ra(errortype, errormsg1, array1, errormsg2, array2, errormsg3)
  integer,            intent(in) :: errortype
  real, dimension(:), intent(in) :: array1
  real, dimension(:), intent(in) :: array2
  character(len=*),      intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3
  character(len=512) :: string

  string = errormsg1//trim(array_to_char(array1))
  string = trim(string)//errormsg2//trim(array_to_char(array2))
  if(present(errormsg3)) string = trim(string)//errormsg3
  call mpp_error_basic( errortype, trim(string))

end subroutine mpp_error_ra_ra
# 254 "mpp/include/mpp_util.inc" 2



!#####################################################################




# 1 "mpp/include/mpp_error_a_s.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_error_ia_is(errortype, errormsg1, array, errormsg2, scalar, errormsg3)
  integer,            intent(in) :: errortype
  integer, dimension(:), intent(in) :: array
  integer,               intent(in) :: scalar
  character(len=*),   intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3

  call mpp_error( errortype, errormsg1, array, errormsg2, (/scalar/), errormsg3)

end subroutine mpp_error_ia_is
# 262 "mpp/include/mpp_util.inc" 2



!#####################################################################




# 1 "mpp/include/mpp_error_a_s.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_error_ia_rs(errortype, errormsg1, array, errormsg2, scalar, errormsg3)
  integer,            intent(in) :: errortype
  integer, dimension(:), intent(in) :: array
  real,               intent(in) :: scalar
  character(len=*),   intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3

  call mpp_error( errortype, errormsg1, array, errormsg2, (/scalar/), errormsg3)

end subroutine mpp_error_ia_rs
# 270 "mpp/include/mpp_util.inc" 2



!#####################################################################




# 1 "mpp/include/mpp_error_a_s.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_error_ra_is(errortype, errormsg1, array, errormsg2, scalar, errormsg3)
  integer,            intent(in) :: errortype
  real, dimension(:), intent(in) :: array
  integer,               intent(in) :: scalar
  character(len=*),   intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3

  call mpp_error( errortype, errormsg1, array, errormsg2, (/scalar/), errormsg3)

end subroutine mpp_error_ra_is
# 278 "mpp/include/mpp_util.inc" 2



!#####################################################################




# 1 "mpp/include/mpp_error_a_s.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_error_ra_rs(errortype, errormsg1, array, errormsg2, scalar, errormsg3)
  integer,            intent(in) :: errortype
  real, dimension(:), intent(in) :: array
  real,               intent(in) :: scalar
  character(len=*),   intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3

  call mpp_error( errortype, errormsg1, array, errormsg2, (/scalar/), errormsg3)

end subroutine mpp_error_ra_rs
# 286 "mpp/include/mpp_util.inc" 2



!#####################################################################




# 1 "mpp/include/mpp_error_s_a.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_error_is_ia(errortype, errormsg1, scalar2, errormsg2, array2, errormsg3)
  integer,            intent(in) :: errortype
  integer,               intent(in) :: scalar2
  integer, dimension(:), intent(in) :: array2
  character(len=*),   intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3

  call mpp_error( errortype, errormsg1, (/scalar2/), errormsg2, array2, errormsg3)

end subroutine mpp_error_is_ia
# 294 "mpp/include/mpp_util.inc" 2



!#####################################################################




# 1 "mpp/include/mpp_error_s_a.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_error_is_ra(errortype, errormsg1, scalar2, errormsg2, array2, errormsg3)
  integer,            intent(in) :: errortype
  integer,               intent(in) :: scalar2
  real, dimension(:), intent(in) :: array2
  character(len=*),   intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3

  call mpp_error( errortype, errormsg1, (/scalar2/), errormsg2, array2, errormsg3)

end subroutine mpp_error_is_ra
# 302 "mpp/include/mpp_util.inc" 2



!#####################################################################




# 1 "mpp/include/mpp_error_s_a.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_error_rs_ia(errortype, errormsg1, scalar2, errormsg2, array2, errormsg3)
  integer,            intent(in) :: errortype
  real,               intent(in) :: scalar2
  integer, dimension(:), intent(in) :: array2
  character(len=*),   intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3

  call mpp_error( errortype, errormsg1, (/scalar2/), errormsg2, array2, errormsg3)

end subroutine mpp_error_rs_ia
# 310 "mpp/include/mpp_util.inc" 2



!#####################################################################




# 1 "mpp/include/mpp_error_s_a.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_error_rs_ra(errortype, errormsg1, scalar2, errormsg2, array2, errormsg3)
  integer,            intent(in) :: errortype
  real,               intent(in) :: scalar2
  real, dimension(:), intent(in) :: array2
  character(len=*),   intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3

  call mpp_error( errortype, errormsg1, (/scalar2/), errormsg2, array2, errormsg3)

end subroutine mpp_error_rs_ra
# 318 "mpp/include/mpp_util.inc" 2



!#####################################################################




# 1 "mpp/include/mpp_error_s_s.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_error_is_is(errortype, errormsg1, scalar1, errormsg2, scalar2, errormsg3)
  integer,            intent(in) :: errortype
  integer, intent(in) :: scalar1
  integer, intent(in) :: scalar2
  character(len=*),   intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3

  call mpp_error( errortype, errormsg1, (/scalar1/), errormsg2, (/scalar2/), errormsg3)

end subroutine mpp_error_is_is
# 326 "mpp/include/mpp_util.inc" 2



!#####################################################################




# 1 "mpp/include/mpp_error_s_s.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_error_is_rs(errortype, errormsg1, scalar1, errormsg2, scalar2, errormsg3)
  integer,            intent(in) :: errortype
  integer, intent(in) :: scalar1
  real, intent(in) :: scalar2
  character(len=*),   intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3

  call mpp_error( errortype, errormsg1, (/scalar1/), errormsg2, (/scalar2/), errormsg3)

end subroutine mpp_error_is_rs
# 334 "mpp/include/mpp_util.inc" 2



!#####################################################################




# 1 "mpp/include/mpp_error_s_s.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_error_rs_is(errortype, errormsg1, scalar1, errormsg2, scalar2, errormsg3)
  integer,            intent(in) :: errortype
  real, intent(in) :: scalar1
  integer, intent(in) :: scalar2
  character(len=*),   intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3

  call mpp_error( errortype, errormsg1, (/scalar1/), errormsg2, (/scalar2/), errormsg3)

end subroutine mpp_error_rs_is
# 342 "mpp/include/mpp_util.inc" 2



!#####################################################################




# 1 "mpp/include/mpp_error_s_s.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_error_rs_rs(errortype, errormsg1, scalar1, errormsg2, scalar2, errormsg3)
  integer,            intent(in) :: errortype
  real, intent(in) :: scalar1
  real, intent(in) :: scalar2
  character(len=*),   intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3

  call mpp_error( errortype, errormsg1, (/scalar1/), errormsg2, (/scalar2/), errormsg3)

end subroutine mpp_error_rs_rs
# 350 "mpp/include/mpp_util.inc" 2



!#####################################################################
function iarray_to_char(iarray) result(string)
integer, intent(in) :: iarray(:)
character(len=256) :: string
character(len=32)  :: chtmp
integer :: i, len_tmp, len_string

 string = ''
 do i=1,size(iarray)
   write(chtmp,'(i16)') iarray(i)
   chtmp = adjustl(chtmp)
   len_tmp = len_trim(chtmp)
   len_string  = len_trim(string)
   string(len_string+1:len_string+len_tmp) = trim(chtmp)
   string(len_string+len_tmp+1:len_string+len_tmp+1) = ','
 enddo
 len_string = len_trim(string)
 string(len_string:len_string) = ' ' ! remove trailing comma

end function iarray_to_char
!#####################################################################
function rarray_to_char(rarray) result(string)
real, intent(in) :: rarray(:)
character(len=256) :: string
character(len=32)  :: chtmp
integer :: i, len_tmp, len_string

 string = ''
 do i=1,size(rarray)
   write(chtmp,'(G16.9)') rarray(i)
   chtmp = adjustl(chtmp)
   len_tmp = len_trim(chtmp)
   len_string  = len_trim(string)
   string(len_string+1:len_string+len_tmp) = trim(chtmp)
   string(len_string+len_tmp+1:len_string+len_tmp+1) = ','
 enddo
 len_string = len_trim(string)
 string(len_string:len_string) = ' ' ! remove trailing comma

end function rarray_to_char

  !> @brief Returns processor ID.
  !!
  !> This returns the unique ID associated with a PE. This number runs
  !! between 0 and <TT>npes-1</TT>, where <TT>npes</TT> is the total
  !! processor count, returned by <TT>mpp_npes</TT>. For a uniprocessor
  !! application this will always return 0.
  function mpp_pe()
    integer :: mpp_pe

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_PE: You must first call mpp_init.' )
    mpp_pe = pe
    return
  end function mpp_pe

  !#####################################################################

  !> @brief Returns processor count for current pelist
  !!
  !> This returns the number of PEs in the current pelist. For a uniprocessor application,
  !! it will always return 1.
  function mpp_npes()
    integer :: mpp_npes

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_NPES: You must first call mpp_init.' )
    mpp_npes = size(peset(current_peset_num)%list(:))
    return
  end function mpp_npes

  !#####################################################################
  function mpp_root_pe()
    integer :: mpp_root_pe

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_ROOT_PE: You must first call mpp_init.' )
    mpp_root_pe = root_pe
    return
  end function mpp_root_pe

  function mpp_comm()
    type(mpi_comm) :: mpp_comm

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'mpp_comm: You must first call mpp_init.' )
    mpp_comm = peset(current_peset_num)%comm
  end function mpp_comm

  function mpp_commID()
    integer :: mpp_commID

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'mpp_commID: You must first call mpp_init.' )
    call mpp_error(NOTE, "mpp_commID() is deprecated. Please use mpp_comm() instead.")
    mpp_commID = peset(current_peset_num)%comm%mpi_val
  end function mpp_commID

  !#####################################################################
  subroutine mpp_set_root_pe(num)
    integer, intent(in) :: num

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SET_ROOT_PE: You must first call mpp_init.' )
    if( .NOT.(ANY(num.EQ.peset(current_peset_num)%list(:))) ) &
         call mpp_error( FATAL, 'MPP_SET_ROOT_PE: you cannot set a root PE outside the current pelist.' )
    root_pe = num
    return
  end subroutine mpp_set_root_pe

  !> @brief Declare a pelist.
  !!
  !> This call is written specifically to accommodate a MPI restriction
  !! that requires a parent communicator to create a child communicator, In
  !! other words: a pelist cannot go off and declare a communicator, but
  !! every PE in the parent, including those not in pelist(:), must get
  !! together for the <TT>MPI_COMM_CREATE</TT> call. The parent is
  !! typically <TT>MPI_COMM_WORLD</TT>, though it could also be a subset
  !! that includes all PEs in <TT>pelist</TT>.
  !!
  !! This call implies synchronization across the PEs in the current
  !! pelist, of which <TT>pelist</TT> is a subset.
  subroutine mpp_declare_pelist_f08( pelist, name, comm )
    integer,                    intent(in) :: pelist(:) !> pelist you are declaring and storing within FMS
    character(len=*), intent(in), optional :: name      !> unique name for an input pelist
    type(mpi_comm), intent(out), optional  :: comm      !> mpi_comm communicator handle
    integer :: i

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_DECLARE_PELIST: You must first call mpp_init.' )
    i = get_peset(pelist)
    write( peset(i)%name,'(a,i2.2)' ) 'PElist', i !default name
    if( PRESENT(name) ) peset(i)%name = name
    if( PRESENT(comm) ) then
      comm = peset(i)%comm
    endif
  end subroutine mpp_declare_pelist_f08

  subroutine mpp_declare_pelist_legacy( pelist, name, commID )
    integer,                    intent(in) :: pelist(:) !> pelist you are declaring and storing within FMS
    character(len=*), intent(in), optional :: name      !> unique name for an input pelist
    integer, intent(out)                   :: commID    !> integral MPI communicator handle
    type(mpi_comm) :: comm

    call mpp_declare_pelist_f08(pelist, name, comm)
    commID = comm%mpi_val
  end subroutine mpp_declare_pelist_legacy

  !#####################################################################

  !> @brief Set context pelist
  !!
  !! This call sets the value of the current pelist, which is the
  !! context for all subsequent "global" calls where the optional
  !! <TT>pelist</TT> argument is omitted. All the PEs that are to be in the
  !! current pelist must call it.
  !!
  !! In MPI, this call may hang unless <TT>pelist</TT> has been previous
  !! declared using @ref mpp_declare_pelist
  !!
  !! If the argument <TT>pelist</TT> is absent, the current pelist is
  !! set to the "world" pelist, of all PEs in the job.
  subroutine mpp_set_current_pelist( pelist, no_sync )
    !Once we branch off into a PE subset, we want subsequent "global" calls to
    !sync only across this subset. This is declared as the current pelist (peset(current_peset_num)%list)
    !when current_peset all pelist ops with no pelist should apply the current pelist.
    !also, we set the start PE in this pelist to be the root_pe.
    !unlike mpp_declare_pelist, this is called by the PEs in the pelist only
    !so if the PEset has not been previously declared, this will hang in MPI.
    !if pelist is omitted, we reset pelist to the world pelist.
    integer, intent(in), optional :: pelist(:)
    logical, intent(in), optional :: no_sync

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SET_CURRENT_PELIST: You must first call mpp_init.' )
    if( PRESENT(pelist) )then
       if( .NOT.ANY(pe.EQ.pelist) )call mpp_error( FATAL, 'MPP_SET_CURRENT_PELIST: pe must be in pelist.' )
       current_peset_num = get_peset(pelist)
    else
       current_peset_num = world_peset_num
    end if
    call mpp_set_root_pe( MINVAL(peset(current_peset_num)%list(:)) )
    if(.not.PRESENT(no_sync))call mpp_sync()  !this is called to make sure everyone in the current pelist is here.
    !      npes = mpp_npes()
    return
  end subroutine mpp_set_current_pelist

  !#####################################################################
  function mpp_get_current_pelist_name()
   ! Simply return the current pelist name
   character(len=len(peset(current_peset_num)%name)) :: mpp_get_current_pelist_name

   mpp_get_current_pelist_name = peset(current_peset_num)%name
  end function mpp_get_current_pelist_name

  !this is created for use by mpp_define_domains within a pelist
  !will be published but not publicized
  subroutine mpp_get_current_pelist_f08( pelist, name, comm )
    integer, intent(out)                    :: pelist(:) !> Array to copy the pelist into
    character(len=*), intent(out), optional :: name      !> Name of the pelist
    type(mpi_comm), intent(out), optional   :: comm      !> mpi_comm communicator handle

    if( size(pelist(:)).NE.size(peset(current_peset_num)%list(:)) ) &
         call mpp_error( FATAL, 'MPP_GET_CURRENT_PELIST: size(pelist) is wrong.' )
    pelist(:) = peset(current_peset_num)%list(:)
    if( PRESENT(name) ) name = peset(current_peset_num)%name
    if( PRESENT(comm) ) then
      comm = peset(current_peset_num)%comm
    endif
  end subroutine mpp_get_current_pelist_f08

  subroutine mpp_get_current_pelist_legacy( pelist, name, commID )
    integer, intent(out)                    :: pelist(:) !> Array to copy the pelist into
    character(len=*), intent(out), optional :: name      !> Name of the pelist
    integer, intent(out)                    :: commID    !> Integral MPI communicator handle
    type(mpi_comm) :: comm

    call mpp_get_current_pelist_f08(pelist, name, comm)
    commID = comm%mpi_val
  end subroutine mpp_get_current_pelist_legacy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                             !
  !                        PERFORMANCE PROFILING CALLS                          !
  !                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief Set the level of granularity of timing measurements.
  !!
  !> This routine and three other routines, mpp_clock_id, mpp_clock_begin(id),
  !!    and mpp_clock_end(id) may be used to time parallel code sections, and
  !!    extract parallel statistics. Clocks are identified by names, which
  !!    should be unique in the first 32 characters. The <TT>mpp_clock_id</TT>
  !!    call initializes a clock of a given name and returns an integer
  !!    <TT>id</TT>. This <TT>id</TT> can be used by subsequent
  !!    <TT>mpp_clock_begin</TT> and <TT>mpp_clock_end</TT> calls set around a
  !!    code section to be timed. Example:
  !!    <PRE>
  !!    integer :: id
  !!    id = mpp_clock_id( 'Atmosphere' )
  !!    call mpp_clock_begin(id)
  !!    call atmos_model()
  !!    call mpp_clock_end()
  !!    </PRE>
  !!     Two flags may be used to alter the behaviour of
  !!     <TT>mpp_clock</TT>. If the flag <TT>MPP_CLOCK_SYNC</TT> is turned on
  !!     by <TT>mpp_clock_id</TT>, the clock calls <TT>mpp_sync</TT> across all
  !!     the PEs in the current pelist at the top of the timed code section,
  !!     but allows each PE to complete the code section (and reach
  !!     <TT>mpp_clock_end</TT>) at different times. This allows us to measure
  !!     load imbalance for a given code section. Statistics are written to
  !!     <TT>stdout</TT> by <TT>mpp_exit</TT>.
  !!
  !!     The flag <TT>MPP_CLOCK_DETAILED</TT> may be turned on by
  !!     <TT>mpp_clock_id</TT> to get detailed communication
  !!     profiles. Communication events of the types <TT>SEND, RECV, BROADCAST,
  !!     REDUCE</TT> and <TT>WAIT</TT> are separately measured for data volume
  !!     and time. Statistics are written to <TT>stdout</TT> by
  !!     <TT>mpp_exit</TT>, and individual PE info is also written to the file
  !!     <TT>mpp_clock.out.####</TT> where <TT>####</TT> is the PE id given by
  !!     <TT>mpp_pe</TT>.
  !!
  !!     The flags <TT>MPP_CLOCK_SYNC</TT> and <TT>MPP_CLOCK_DETAILED</TT> are
  !!     integer parameters available by use association, and may be summed to
  !!     turn them both on.
  !!
  !!     While the nesting of clocks is allowed, please note that turning on
  !!     the non-optional flags on inner clocks has certain subtle issues.
  !!     Turning on <TT>MPP_CLOCK_SYNC</TT> on an inner
  !!     clock may distort outer clock measurements of load imbalance. Turning
  !!     on <TT>MPP_CLOCK_DETAILED</TT> will stop detailed measurements on its
  !!     outer clock, since only one detailed clock may be active at one time.
  !!     Also, detailed clocks only time a certain number of events per clock
  !!     (currently 40000) to conserve memory. If this array overflows, a
  !!     warning message is printed, and subsequent events for this clock are
  !!     not timed.
  !!
  !!     Timings are done using the <TT>f90</TT> standard
  !!     <TT>system_clock_default</TT> intrinsic.
  !!
  !!     The resolution of system_clock_default is often too coarse for use except
  !!     across large swaths of code. On SGI systems this is transparently
  !!     overloaded with a higher resolution clock made available in a
  !!     non-portable fortran interface made available by
  !!     <TT>nsclock.c</TT>. This approach will eventually be extended to other
  !!     platforms.
  !!
  !!     New behaviour added at the Havana release allows the user to embed
  !!     profiling calls at varying levels of granularity all over the code,
  !!     and for any particular run, set a threshold of granularity so that
  !!     finer-grained clocks become dormant.
  !!
  !!     The threshold granularity is held in the private module variable
  !!     <TT>clock_grain</TT>. This value may be modified by the call
  !!     <TT>mpp_clock_set_grain</TT>, and affect clocks initiated by
  !!     subsequent calls to <TT>mpp_clock_id</TT>. The value of
  !!     <TT>clock_grain</TT> is set to an arbitrarily large number initially.
  !!
  !!     Clocks initialized by <TT>mpp_clock_id</TT> can set a new optional
  !!     argument <TT>grain</TT> setting their granularity level. Clocks check
  !!     this level against the current value of <TT>clock_grain</TT>, and are
  !!     only triggered if they are <I>at or below ("coarser than")</I> the
  !!     threshold. Finer-grained clocks are dormant for that run.
  !!
  !!The following grain levels are pre-defined:
  !!
  !!<pre>
  !!!predefined clock granularities, but you can use any integer
  !!!using CLOCK_LOOP and above may distort coarser-grain measurements
  !!  integer, parameter, public :: CLOCK_COMPONENT=1 !component level, e.g model, exchange
  !!  integer, parameter, public :: CLOCK_SUBCOMPONENT=11 !top level within a model component, e.g dynamics, physics
  !!  integer, parameter, public :: CLOCK_MODULE=21 !module level, e.g main subroutine of a physics module
  !!  integer, parameter, public :: CLOCK_ROUTINE=31 !level of individual subroutine or function
  !!  integer, parameter, public :: CLOCK_LOOP=41 !loops or blocks within a routine
  !!  integer, parameter, public :: CLOCK_INFRA=51 !infrastructure level, e.g halo update
  !!</pre>
  !!
  !!     Note that subsequent changes to <TT>clock_grain</TT> do not
  !!     change the status of already initiated clocks, and that if the
  !!     optional <TT>grain</TT> argument is absent, the clock is always
  !!     triggered. This guarantees backward compatibility.
  subroutine mpp_clock_set_grain( grain )
    integer, intent(in) :: grain
    !set the granularity of times: only clocks whose grain is lower than
    !clock_grain are triggered, finer-grained clocks are dormant.
    !clock_grain is initialized to CLOCK_LOOP, so all clocks above the loop level
    !are triggered if this is never called.
    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_CLOCK_SET_GRAIN: You must first call mpp_init.' )

    clock_grain = grain
    return
  end subroutine mpp_clock_set_grain

  !#####################################################################
  subroutine clock_init( id, name, flags, grain )
    integer,           intent(in) :: id
    character(len=*),  intent(in) :: name
    integer, intent(in), optional :: flags, grain
    integer                       :: i

    clocks(id)%name = name
    clocks(id)%hits = 0
    clocks(id)%tick = 0
    clocks(id)%total_ticks = 0
    clocks(id)%sync_on_begin = .FALSE.
    clocks(id)%detailed      = .FALSE.
    clocks(id)%peset_num = current_peset_num
    if( PRESENT(flags) )then
       if( BTEST(flags,0) )clocks(id)%sync_on_begin = .TRUE.
       if( BTEST(flags,1) )clocks(id)%detailed      = .TRUE.
    end if
    clocks(id)%grain = 0
    if( PRESENT(grain) )clocks(id)%grain = grain
    if( clocks(id)%detailed )then
       allocate( clocks(id)%events(MAX_EVENT_TYPES) )
       clocks(id)%events(EVENT_ALLREDUCE)%name = 'ALLREDUCE'
       clocks(id)%events(EVENT_BROADCAST)%name = 'BROADCAST'
       clocks(id)%events(EVENT_RECV)%name = 'RECV'
       clocks(id)%events(EVENT_SEND)%name = 'SEND'
       clocks(id)%events(EVENT_WAIT)%name = 'WAIT'
       do i=1,MAX_EVENT_TYPES
          clocks(id)%events(i)%ticks(:) = 0
          clocks(id)%events(i)%bytes(:) = 0
          clocks(id)%events(i)%calls = 0
       end do
       clock_summary(id)%name = name
       clock_summary(id)%event(EVENT_ALLREDUCE)%name = 'ALLREDUCE'
       clock_summary(id)%event(EVENT_BROADCAST)%name = 'BROADCAST'
       clock_summary(id)%event(EVENT_RECV)%name = 'RECV'
       clock_summary(id)%event(EVENT_SEND)%name = 'SEND'
       clock_summary(id)%event(EVENT_WAIT)%name = 'WAIT'
       do i=1,MAX_EVENT_TYPES
          clock_summary(id)%event(i)%msg_size_sums(:) = 0.0
          clock_summary(id)%event(i)%msg_time_sums(:) = 0.0
          clock_summary(id)%event(i)%total_data = 0.0
          clock_summary(id)%event(i)%total_time = 0.0
          clock_summary(id)%event(i)%msg_size_cnts(:) = 0
          clock_summary(id)%event(i)%total_cnts = 0
       end do
    end if
    return
  end subroutine clock_init

  !#####################################################################
  !> Return an ID for a new or existing clock
  function mpp_clock_id( name, flags, grain )
    integer                       :: mpp_clock_id
    character(len=*),  intent(in) :: name
    integer, intent(in), optional :: flags, grain

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_CLOCK_ID: You must first call mpp_init.')

    !if grain is present, the clock is only triggered if it
    !is low ("coarse") enough: compared to clock_grain
    !finer-grained clocks are dormant.
    !if grain is absent, clock is triggered.
    if( PRESENT(grain) )then
       if( grain.GT.clock_grain )then
          mpp_clock_id = 0
          return
       end if
    end if
    mpp_clock_id = 1

    if( clock_num.EQ.0 )then  !first
       clock_num = mpp_clock_id
       call clock_init(mpp_clock_id,name,flags)
    else
       FIND_CLOCK: do while( trim(name).NE.trim(clocks(mpp_clock_id)%name) )
          mpp_clock_id = mpp_clock_id + 1
          if( mpp_clock_id.GT.clock_num )then
             if( mpp_clock_id.GT.MAX_CLOCKS )then
                call mpp_error( FATAL, 'MPP_CLOCK_ID: too many clock requests, ' // &
                      'check your clock id request or increase MAX_CLOCKS.')
             else               !new clock: initialize
                clock_num = mpp_clock_id
                call clock_init(mpp_clock_id,name,flags,grain)
                exit FIND_CLOCK
             end if
          end if
       end do FIND_CLOCK
    endif
    return
  end function mpp_clock_id

  !#####################################################################
  subroutine mpp_clock_begin(id)
    integer, intent(in) :: id

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_CLOCK_BEGIN: You must first call mpp_init.' )
    if( .not. mpp_record_timing_data)return
    if( id.EQ.0 )return
    if( id.LT.0 .OR. id.GT.clock_num )call mpp_error( FATAL, 'MPP_CLOCK_BEGIN: invalid id.' )

!$OMP MASTER
    if( clocks(id)%peset_num.NE.current_peset_num ) &
         call mpp_error( FATAL, 'MPP_CLOCK_BEGIN: cannot change pelist context of a clock.' )
    if( clocks(id)%is_on) call mpp_error(FATAL, 'MPP_CLOCK_BEGIN: mpp_clock_begin is called again '// &
                'before calling mpp_clock_end for the clock '//trim(clocks(id)%name) )
    if( clocks(id)%sync_on_begin .OR. sync_all_clocks )then
       !do an untimed sync at the beginning of the clock
       !this puts all PEs in the current pelist on par, so that measurements begin together
       !ending time will be different, thus measuring load imbalance for this clock.
       call mpp_sync()
    end if

    if (debug) then
      num_clock_ids = num_clock_ids+1
      if(num_clock_ids > MAX_CLOCKS)call mpp_error(FATAL,'MPP_CLOCK_BEGIN: max num previous_clock exceeded.' )
      previous_clock(num_clock_ids) = current_clock
      current_clock = id
    endif
    call system_clock_default( clocks(id)%tick )
    clocks(id)%hits = clocks(id)%hits + 1
    clocks(id)%is_on = .true.
!$OMP END MASTER
    return
  end subroutine mpp_clock_begin

  !#####################################################################
  subroutine mpp_clock_end(id)
    integer, intent(in) :: id
    integer(i8_kind)  :: delta
    integer             :: errunit

    if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_CLOCK_END: You must first call mpp_init.' )
    if( .not. mpp_record_timing_data)return
    if( id.EQ.0 )return
    if( id.LT.0 .OR. id.GT.clock_num )call mpp_error( FATAL, 'MPP_CLOCK_BEGIN: invalid id.' )
!$OMP MASTER
    if( .NOT. clocks(id)%is_on) call mpp_error(FATAL, 'MPP_CLOCK_END: mpp_clock_end is called '// &
                'before calling mpp_clock_begin for the clock '//trim(clocks(id)%name) )

    call system_clock_default(end_tick)
    if( clocks(id)%peset_num.NE.current_peset_num ) &
         call mpp_error( FATAL, 'MPP_CLOCK_END: cannot change pelist context of a clock.' )
    delta = end_tick - clocks(id)%tick
    if( delta.LT.0 )then
       errunit = stderr()
       write( errunit,* )'pe, id, start_tick, end_tick, delta, max_ticks=', pe, id, clocks(id)%tick, end_tick, &
            &  delta, max_ticks
       delta = delta + max_ticks + 1
       call mpp_error( WARNING, 'MPP_CLOCK_END: Clock rollover, assumed single roll.' )
    end if
    clocks(id)%total_ticks = clocks(id)%total_ticks + delta
    if (debug) then
      if(num_clock_ids < 1) call mpp_error(NOTE,'MPP_CLOCK_END: min num previous_clock < 1.' )
      current_clock = previous_clock(num_clock_ids)
      num_clock_ids = num_clock_ids-1
    endif
    clocks(id)%is_on = .false.
!$OMP END MASTER
    return
  end subroutine mpp_clock_end

 !#####################################################################
  subroutine mpp_record_time_start()

     mpp_record_timing_data = .TRUE.

  end subroutine mpp_record_time_start

  !#####################################################################
  subroutine mpp_record_time_end()

     mpp_record_timing_data = .FALSE.

  end subroutine mpp_record_time_end


  !#####################################################################
  subroutine increment_current_clock( event_id, bytes )
    integer,           intent(in) :: event_id
    integer, intent(in), optional :: bytes
    integer                       :: n
    integer(i8_kind)            :: delta
    integer                       :: errunit

    if( .not. mpp_record_timing_data )return
    if( .not.debug .or. (current_clock.EQ.0) )return
    if( current_clock.LT.0 .OR. current_clock.GT.clock_num )call mpp_error( FATAL, &
       &  'MPP_CLOCK_BEGIN: invalid current_clock.' )
    if( .NOT.clocks(current_clock)%detailed )return
    call system_clock_default(end_tick)
    n = clocks(current_clock)%events(event_id)%calls + 1

    if( n.EQ.MAX_EVENTS )call mpp_error( WARNING, &
         'MPP_CLOCK: events exceed MAX_EVENTS, ignore detailed profiling data for clock '// &
         & trim(clocks(current_clock)%name) )
    if( n.GT.MAX_EVENTS )return

    clocks(current_clock)%events(event_id)%calls = n
    delta = end_tick - start_tick
    if( delta.LT.0 )then
       errunit = stderr()
       write( errunit,* )'pe, event_id, start_tick, end_tick, delta, max_ticks=', &
                           pe, event_id, start_tick, end_tick, delta, max_ticks
       delta = delta + max_ticks + 1
       call mpp_error( WARNING, 'MPP_CLOCK_END: Clock rollover, assumed single roll.' )
    end if
    clocks(current_clock)%events(event_id)%ticks(n) = delta
    if( PRESENT(bytes) )clocks(current_clock)%events(event_id)%bytes(n) = bytes
    return
  end subroutine increment_current_clock

  !#####################################################################

  subroutine dump_clock_summary()

    real              :: total_time,total_time_all,total_data
    real              :: msg_size,eff_BW,s
    integer           :: SD_UNIT, total_calls
    integer           :: j,k,ct, msg_cnt
    character(len=2)  :: u
    character(len=FMS_FILE_LEN) :: filename
    character(len=20),dimension(MAX_BINS),save :: bin

    data bin( 1)  /'  0   -    8    B:  '/
    data bin( 2)  /'  8   -   16    B:  '/
    data bin( 3)  /' 16   -   32    B:  '/
    data bin( 4)  /' 32   -   64    B:  '/
    data bin( 5)  /' 64   -  128    B:  '/
    data bin( 6)  /'128   -  256    B:  '/
    data bin( 7)  /'256   -  512    B:  '/
    data bin( 8)  /'512   - 1024    B:  '/
    data bin( 9)  /'  1.0 -    2.1 KB:  '/
    data bin(10)  /'  2.1 -    4.1 KB:  '/
    data bin(11)  /'  4.1 -    8.2 KB:  '/
    data bin(12)  /'  8.2 -   16.4 KB:  '/
    data bin(13)  /' 16.4 -   32.8 KB:  '/
    data bin(14)  /' 32.8 -   65.5 KB:  '/
    data bin(15)  /' 65.5 -  131.1 KB:  '/
    data bin(16)  /'131.1 -  262.1 KB:  '/
    data bin(17)  /'262.1 -  524.3 KB:  '/
    data bin(18)  /'524.3 - 1048.6 KB:  '/
    data bin(19)  /'  1.0 -    2.1 MB:  '/
    data bin(20)  /' >2.1          MB:  '/

    if( .NOT.ANY(clocks(1:clock_num)%detailed) )return
    write( filename,'(a,i6.6)' )'mpp_clock.out.', pe

    open(newunit=SD_UNIT,file=trim(filename),form='formatted')

    COMM_TYPE: do ct = 1,clock_num

       if( .NOT.clocks(ct)%detailed )cycle
       write(SD_UNIT,*) &
            clock_summary(ct)%name(1:15),' Communication Data for PE ',pe

       write(SD_UNIT,*) ' '
       write(SD_UNIT,*) ' '

       total_time_all = 0.0
       EVENT_TYPE: do k = 1,MAX_EVENT_TYPES-1

          if(clock_summary(ct)%event(k)%total_time == 0.0)cycle

          total_time = clock_summary(ct)%event(k)%total_time
          total_time_all = total_time_all + total_time
          total_data = clock_summary(ct)%event(k)%total_data
          total_calls = int(clock_summary(ct)%event(k)%total_cnts)

          write(SD_UNIT,1000) clock_summary(ct)%event(k)%name(1:9) // ':'

          write(SD_UNIT,1001) 'Total Data: ',total_data*1.0e-6, &
               'MB; Total Time: ', total_time, &
               'secs; Total Calls: ',total_calls

          write(SD_UNIT,*) ' '
          write(SD_UNIT,1002) '     Bin            Counts      Avg Size        Eff B/W'
          write(SD_UNIT,*) ' '

          BIN_LOOP: do j=1,MAX_BINS

             if(clock_summary(ct)%event(k)%msg_size_cnts(j)==0)cycle

             if(j<=8)then
                s = 1.0
                u = ' B'
             elseif(j<=18)then
                s = 1.0e-3
                u = 'KB'
             else
                s = 1.0e-6
                u = 'MB'
             endif

             msg_cnt = int(clock_summary(ct)%event(k)%msg_size_cnts(j))
             msg_size = &
                  s*(clock_summary(ct)%event(k)%msg_size_sums(j)/real(msg_cnt))
             eff_BW = (1.0e-6)*( clock_summary(ct)%event(k)%msg_size_sums(j) / &
                  clock_summary(ct)%event(k)%msg_time_sums(j) )

             write(SD_UNIT,1003) bin(j),msg_cnt,msg_size,u,eff_BW

          end do BIN_LOOP

          write(SD_UNIT,*) ' '
          write(SD_UNIT,*) ' '
       end do EVENT_TYPE

       ! "Data-less" WAIT

       if(clock_summary(ct)%event(MAX_EVENT_TYPES)%total_time>0.0)then

          total_time = clock_summary(ct)%event(MAX_EVENT_TYPES)%total_time
          total_time_all = total_time_all + total_time
          total_calls = int(clock_summary(ct)%event(MAX_EVENT_TYPES)%total_cnts)

          write(SD_UNIT,1000) clock_summary(ct)%event(MAX_EVENT_TYPES)%name(1:9) // ':'

          write(SD_UNIT,1004) 'Total Calls: ',total_calls,'; Total Time: ', &
               total_time,'secs'

       endif

       write(SD_UNIT,*) ' '
       write(SD_UNIT,1005) 'Total communication time spent for ' // &
            clock_summary(ct)%name(1:9) // ': ',total_time_all,'secs'
       write(SD_UNIT,*) ' '
       write(SD_UNIT,*) ' '
       write(SD_UNIT,*) ' '

    end do COMM_TYPE

    close(SD_UNIT)

1000 format(a)
1001 format(a,f8.2,a,f8.2,a,i6)
1002 format(a)
1003 format(a,i6,'    ','  ',f9.1,a,'    ',f9.2,'MB/sec')
1004 format(a,i8,a,f9.2,a)
1005 format(a,f9.2,a)
    return
  end subroutine dump_clock_summary

  !#####################################################################

  integer function get_unit()

    integer,save :: i
    logical      :: l_open

    if (pe == root_pe) call mpp_error(WARNING, &
        'get_unit is deprecated and will be removed in a future release, please use the Fortran intrinsic newunit')
    do i=10,99
       inquire(unit=i,opened=l_open)
       if(.not.l_open)exit
    end do

    if(i==100)then
       call mpp_error(FATAL,'Unable to get I/O unit')
    else
       get_unit = i
    endif

    return
  end function get_unit

  !#####################################################################

  subroutine sum_clock_data()

    integer :: i,j,k,ct,event_size,event_cnt
    real    :: msg_time

    CLOCK_TYPE: do ct=1,clock_num
       if( .NOT.clocks(ct)%detailed )cycle
       EVENT_TYPE: do j=1,MAX_EVENT_TYPES-1
          event_cnt = clocks(ct)%events(j)%calls
          EVENT_SUMMARY: do i=1,event_cnt

             clock_summary(ct)%event(j)%total_cnts = &
                  clock_summary(ct)%event(j)%total_cnts + 1

             event_size = int(clocks(ct)%events(j)%bytes(i))

             k = find_bin(event_size)

             clock_summary(ct)%event(j)%msg_size_cnts(k) = &
                  clock_summary(ct)%event(j)%msg_size_cnts(k) + 1

             clock_summary(ct)%event(j)%msg_size_sums(k) = &
                  clock_summary(ct)%event(j)%msg_size_sums(k) &
                  + clocks(ct)%events(j)%bytes(i)

             clock_summary(ct)%event(j)%total_data = &
                  clock_summary(ct)%event(j)%total_data &
                  + clocks(ct)%events(j)%bytes(i)

             msg_time = clocks(ct)%events(j)%ticks(i)
             msg_time = tick_rate * real( clocks(ct)%events(j)%ticks(i) )

             clock_summary(ct)%event(j)%msg_time_sums(k) = &
                  clock_summary(ct)%event(j)%msg_time_sums(k) + msg_time

             clock_summary(ct)%event(j)%total_time = &
                  clock_summary(ct)%event(j)%total_time + msg_time

          end do EVENT_SUMMARY
       end do EVENT_TYPE

       j = MAX_EVENT_TYPES ! WAITs
       ! "msg_size_cnts" doesn't really mean anything for WAIT
       ! but position will be used to store number of counts for now.

       event_cnt = clocks(ct)%events(j)%calls
       clock_summary(ct)%event(j)%msg_size_cnts(1) = event_cnt
       clock_summary(ct)%event(j)%total_cnts       = event_cnt

       msg_time = tick_rate * real( sum ( clocks(ct)%events(j)%ticks(1:event_cnt) ) )
       clock_summary(ct)%event(j)%msg_time_sums(1) = &
            clock_summary(ct)%event(j)%msg_time_sums(1) + msg_time

       clock_summary(ct)%event(j)%total_time = clock_summary(ct)%event(j)%msg_time_sums(1)

    end do CLOCK_TYPE

    return
  contains
    integer function find_bin(event_size)

      integer,intent(in) :: event_size
      integer            :: k,msg_size

      msg_size = 8
      k = 1
      do while(event_size>msg_size .and. k<MAX_BINS)
         k = k+1
         msg_size = msg_size*2
      end do
      find_bin = k
      return
    end function find_bin

  end subroutine sum_clock_data

  !#####################################################################
  !> This routine will double the size of peset and copy the original peset data
  !! into the expanded one. The maximum allowed to expand is PESET_MAX.
  subroutine expand_peset()
     integer :: old_peset_max,n
     type(communicator), allocatable :: peset_old(:)

     old_peset_max = current_peset_max
     if(old_peset_max .GE. PESET_MAX) call mpp_error(FATAL, &
         "mpp_mod(expand_peset): the number of peset reached PESET_MAX, increase PESET_MAX or contact developer")

     ! copy data to a tempoary data
     allocate(peset_old(0:old_peset_max))
     do n = 0, old_peset_max
        peset_old(n)%count      = peset(n)%count
        peset_old(n)%comm       = peset(n)%comm
        peset_old(n)%group      = peset(n)%group
        peset_old(n)%name       = peset(n)%name
        peset_old(n)%start      = peset(n)%start
        peset_old(n)%log2stride = peset(n)%log2stride

        if( ASSOCIATED(peset(n)%list) ) then
           allocate(peset_old(n)%list(size(peset(n)%list(:))) )
           peset_old(n)%list(:) = peset(n)%list(:)
           deallocate(peset(n)%list)
        endif
     enddo
     deallocate(peset)

     ! create the new peset
     current_peset_max = min(PESET_MAX, 2*old_peset_max)
     allocate(peset(0:current_peset_max))
     peset(:)%count = -1
     peset(:)%comm = mpi_comm_null
     peset(:)%group = mpi_group_null
     peset(:)%start = -1
     peset(:)%log2stride = -1
     peset(:)%name  = " "
     do n = 0, old_peset_max
        peset(n)%count      = peset_old(n)%count
        peset(n)%comm       = peset_old(n)%comm
        peset(n)%group      = peset_old(n)%group
        peset(n)%name       = peset_old(n)%name
        peset(n)%start      = peset_old(n)%start
        peset(n)%log2stride = peset_old(n)%log2stride

        if( ASSOCIATED(peset_old(n)%list) ) then
           allocate(peset(n)%list(size(peset_old(n)%list(:))) )
           peset(n)%list(:) = peset_old(n)%list(:)
           deallocate(peset_old(n)%list)
        endif
     enddo
     deallocate(peset_old)

     call mpp_error(NOTE, "mpp_mod(expand_peset): size of peset is expanded to ", current_peset_max)

  end subroutine expand_peset
  !#####################################################################

  function uppercase (cs)
    character(len=*), intent(in) :: cs
    character(len=len(cs)),target       :: uppercase
    integer                      :: k,tlen
    character, pointer :: ca
    integer, parameter :: co=iachar('A')-iachar('a') ! case offset
    !The transfer function truncates the string with xlf90_r
    tlen = len_trim(cs)
    if(tlen <= 0) then      ! catch IBM compiler bug
       uppercase = cs  ! simply return input blank string
    else
    uppercase = cs(1:tlen)
    do k=1, tlen
       ca => uppercase(k:k)
       if(ca >= "a" .and. ca <= "z") ca = achar(ichar(ca)+co)
    enddo
    endif
  end function uppercase

!#######################################################################

  function lowercase (cs)
    character(len=*), intent(in) :: cs
    character(len=len(cs)),target       :: lowercase
    integer, parameter :: co=iachar('a')-iachar('A') ! case offset
    integer                        :: k,tlen
    character, pointer :: ca
!  The transfer function truncates the string with xlf90_r
    tlen = len_trim(cs)
    if(tlen <= 0) then      ! catch IBM compiler bug
       lowercase = cs  ! simply return input blank string
    else
    lowercase = cs(1:tlen)
    do k=1, tlen
       ca => lowercase(k:k)
       if(ca >= "A" .and. ca <= "Z") ca = achar(ichar(ca)+co)
    enddo
    endif
  end function lowercase


  !#######################################################################

!-----------------------------------------------------------------------
!
! AUTHOR: Rusty Benson (rusty.benson@noaa.gov)
!
!
! THESE LINES MUST BE PRESENT IN MPP.F90
!
! ! public variable needed for reading an input nml file from an internal file
!   character(len=:), dimension(:), allocatable, public :: input_nml_file
!

!-----------------------------------------------------------------------

!> Reads an existing input nml file into a character array and broadcasts
!! it to the non-root mpi-tasks. This allows the use of reads from an
!! internal file for namelist settings (requires 2003 compliant compiler)
!!
!! read(input_nml_file, nml=<name_nml>, iostat=status)
!!
!!
  subroutine read_input_nml(pelist_name_in, alt_input_nml_path)

! Include variable "version" to be written to log file.

# 1 "./include/file_version.h" 1
! -*-f90-*-
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************




  character(len=*), parameter :: version = 'unknown'
# 1248 "mpp/include/mpp_util.inc" 2

    character(len=*), intent(in), optional :: pelist_name_in
    character(len=*), intent(in), optional :: alt_input_nml_path
! private variables
    integer :: log_unit
    integer :: i
    integer, dimension(2) :: lines_and_length
    logical :: file_exist
    character(len=len(peset(current_peset_num)%name)) :: pelist_name
    character(len=FMS_PATH_LEN) :: filename

! check the status of input_nml_file
    if ( allocated(input_nml_file) ) then
      deallocate(input_nml_file)
    endif

! the following code is necessary for using alternate namelist files (nests, stretched grids, etc)
    if (PRESENT(pelist_name_in)) then
      ! test to make sure length of pelist_name_in is <= pelist_name
      if (LEN(pelist_name_in) > LEN(pelist_name)) then
        call mpp_error(FATAL,  &
           "mpp_util.inc: read_input_nml optional argument pelist_name_in has size greater than local pelist_name")
      else
        pelist_name = pelist_name_in
      endif
    else
      pelist_name = mpp_get_current_pelist_name()
    endif
    filename='input_'//trim(pelist_name)//'.nml'
    inquire(FILE=filename, EXIST=file_exist)
    if (.not. file_exist ) then
       if (present(alt_input_nml_path)) then
          filename = alt_input_nml_path
       else
          filename = 'input.nml'
       end if
    endif
    lines_and_length = get_ascii_file_num_lines_and_length(filename)
    allocate(character(len=lines_and_length(2))::input_nml_file(lines_and_length(1)))
    call read_ascii_file(filename, lines_and_length(2), input_nml_file)

! write info logfile
    if (pe == root_pe) then
       log_unit = stdlog()
       write(log_unit,'(a)')  '========================================================================'
       write(log_unit,'(a)')  'READ_INPUT_NML: '//trim(version)
       write(log_unit,'(a)')  'READ_INPUT_NML: '//trim(filename)//' '
       do i = 1, lines_and_length(1)
          write(log_unit,*) trim(input_nml_file(i))
       enddo
    end if
  end subroutine read_input_nml


  !#######################################################################
  !z1l: This is extracted from read_ascii_file
  function get_ascii_file_num_lines(FILENAME, LENGTH, PELIST)
    character(len=*), intent(in) :: FILENAME
    integer, intent(in) :: LENGTH
    integer, intent(in), optional, dimension(:) :: PELIST

    integer :: num_lines, get_ascii_file_num_lines
    character(len=LENGTH) :: str_tmp
    character(len=5) :: text
    integer :: status, f_unit, from_pe
    logical :: file_exist

    if( read_ascii_file_on) then
       call mpp_error(FATAL,  &
          "mpp_util.inc: get_ascii_file_num_lines is called again before calling read_ascii_file")
    endif
    read_ascii_file_on = .true.

    from_pe = root_pe
    get_ascii_file_num_lines = -1
    num_lines = -1
    if ( pe == root_pe ) then
       inquire(FILE=FILENAME, EXIST=file_exist)

       if ( file_exist ) then
          open(newunit=f_unit, FILE=FILENAME, ACTION='READ', STATUS='OLD', IOSTAT=status)

          if ( status .ne. 0 ) then
             write (UNIT=text, FMT='(I5)') status
             call mpp_error(FATAL, 'get_ascii_file_num_lines: Error opening file:' //trim(FILENAME)// &
                            '.  (IOSTAT = '//trim(text)//')')
          else
             num_lines = 1
             do
                read (UNIT=f_unit, FMT='(A)', IOSTAT=status) str_tmp
                if ( status .lt. 0 ) then
                   ! deprecate num_lines by 1 and ensure num_lines is at least 1
                   num_lines = max (num_lines - 1, 1)
                   exit
                endif
                if ( status .gt. 0 ) then
                   write (UNIT=text, FMT='(I5)') num_lines
                   call mpp_error(FATAL, 'get_ascii_file_num_lines: Error reading line '//trim(text)// &
                        ' in file '//trim(FILENAME)//'.')
                end if
                if ( len_trim(str_tmp) == LENGTH ) then
                   write(UNIT=text, FMT='(I5)') length
                   call mpp_error(FATAL, 'get_ascii_file_num_lines: Length of output string ('//trim(text)//&
                                       & ' is too small. Increase the LENGTH value.')
                end if
                num_lines = num_lines + 1
             end do
             close(UNIT=f_unit)
          end if
       else
          call mpp_error(FATAL, 'get_ascii_file_num_lines: File '//trim(FILENAME)//' does not exist.')
       end if
    end if

    ! Broadcast number of lines
    call mpp_broadcast(num_lines, from_pe, PELIST=PELIST)
    get_ascii_file_num_lines = num_lines

  end function get_ascii_file_num_lines

  !#######################################################################
  !> @brief Function to determine the maximum line length and number of lines from an ascii file
  function get_ascii_file_num_lines_and_length(FILENAME, PELIST)
    character(len=*), intent(in) :: FILENAME !< name of the file to be read
    integer, intent(in), optional, dimension(:) :: PELIST !< optional pelist

    integer, dimension(2) :: get_ascii_file_num_lines_and_length !< number of lines (1) and
                                                                 !! max line length (2)
    integer :: num_lines, max_length
    integer, parameter :: LENGTH=1024
    character(len=LENGTH) :: str_tmp
    character(len=5) :: text
    integer :: status, f_unit, from_pe
    logical :: file_exist

    if( read_ascii_file_on) then
       call mpp_error(FATAL,  &
          "mpp_util.inc: get_ascii_file_num_lines is called again before calling read_ascii_file")
    endif
    read_ascii_file_on = .true.

    from_pe = root_pe
    get_ascii_file_num_lines_and_length = -1
    num_lines = -1
    max_length = -1
    if ( pe == root_pe ) then
       inquire(FILE=FILENAME, EXIST=file_exist)

       if ( file_exist ) then
          open(newunit=f_unit, FILE=FILENAME, ACTION='READ', STATUS='OLD', IOSTAT=status)

          if ( status .ne. 0 ) then
             write (UNIT=text, FMT='(I5)') status
             call mpp_error(FATAL, 'get_ascii_file_num_lines: Error opening file:' //trim(FILENAME)// &
                            '.  (IOSTAT = '//trim(text)//')')
          else
             num_lines = 1
             max_length = 1
             do
                read (UNIT=f_unit, FMT='(A)', IOSTAT=status) str_tmp
                if ( status .lt. 0 ) then
                   ! deprecate num_lines by 1 and ensure num_lines is at least 1
                   num_lines = max (num_lines - 1, 1)
                   exit
                endif
                if ( status .gt. 0 ) then
                   write (UNIT=text, FMT='(I5)') num_lines
                   call mpp_error(FATAL, 'get_ascii_file_num_lines: Error reading line '//trim(text)// &
                        ' in file '//trim(FILENAME)//'.')
                end if
                if ( len_trim(str_tmp) == LENGTH) then
                   write(UNIT=text, FMT='(I5)') LENGTH
                   call mpp_error(FATAL, 'get_ascii_file_num_lines: Length of output string ('//trim(text)//&
                                       & ' is too small. Increase the LENGTH value.')
                end if
                if (len_trim(str_tmp) > max_length) max_length = len_trim(str_tmp)
                num_lines = num_lines + 1
             end do
             close(UNIT=f_unit)
          end if
       else
          call mpp_error(FATAL, 'get_ascii_file_num_lines: File '//trim(FILENAME)//' does not exist.')
       end if
       max_length = max_length+1
    end if

    ! Broadcast number of lines
    call mpp_broadcast(num_lines, from_pe, PELIST=PELIST)
    call mpp_broadcast(max_length, from_pe, PELIST=PELIST)
    get_ascii_file_num_lines_and_length(1) = num_lines
    get_ascii_file_num_lines_and_length(2) = max_length

  end function get_ascii_file_num_lines_and_length

  !-----------------------------------------------------------------------
  !
  ! AUTHOR: Rusty Benson <rusty.benson@noaa.gov>,
  !         Seth Underwood <Seth.Underwood@noaa.gov>
  !
  !-----------------------------------------------------------------------
  ! subroutine READ_ASCII_FILE
  !
  !
  !> Reads any ascii file into a character array and broadcasts
  !! it to the non-root mpi-tasks.  Based off READ_INPUT_NML.
  !!
  !! Passed in 'Content' array, must be of the form:
  !! character(len=LENGTH), dimension(:), allocatable :: array_name
  !!
  !! Reads from this array must be done in a do loop over the number of
  !! lines, i.e.:
  !!
  !! do i=1, num_lines
  !!    read (UNIT=array_name(i), FMT=*) var1, var2, ...
  !! end do
  subroutine read_ascii_file(FILENAME, LENGTH, Content, PELIST)
    character(len=*),    intent(in)               :: FILENAME
    integer,             intent(in)               :: LENGTH
    character(len=*), intent(inout), dimension(:) :: Content
    integer, intent(in), optional,   dimension(:) :: PELIST

    ! Include variable "version" to be written to log file.

# 1 "./include/file_version.h" 1
! -*-f90-*-
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************




  character(len=*), parameter :: version = 'unknown'
# 1471 "mpp/include/mpp_util.inc" 2

    character(len=5) :: text
    logical :: file_exist
    integer :: status, f_unit, log_unit
    integer :: from_pe
    integer :: pnum_lines, num_lines
    character(len=LENGTH) :: str_tmp !< Temporary variable to store line from file

    if( .NOT. read_ascii_file_on) then
       call mpp_error(FATAL,  &
          "mpp_util.inc: get_ascii_file_num_lines needs to be called before calling read_ascii_file")
    endif
    read_ascii_file_on = .false.

    from_pe = root_pe
    num_lines = size(Content(:))

    if ( pe == root_pe ) then
       ! write info logfile
       log_unit = stdlog()
       write(log_unit,'(a)')  '========================================================================'
       write(log_unit,'(a)')  'READ_ASCII_FILE: '//trim(version)
       write(log_unit,'(a)')  'READ_ASCII_FILE: File: '//trim(FILENAME)

       inquire(FILE=FILENAME, EXIST=file_exist)

       if ( file_exist ) then
          open(newunit=f_unit, FILE=FILENAME, ACTION='READ', STATUS='OLD', IOSTAT=status)

          if ( status .ne. 0 ) then
             write (UNIT=text, FMT='(I5)') status
             call mpp_error(FATAL, 'READ_ASCII_FILE: Error opening file: '// &
                            & trim(FILENAME)//'.  (IOSTAT = '//trim(text)//')')
          else

             if ( num_lines .gt. 0 ) then
                Content(:) = ' '

                rewind(UNIT=f_unit, IOSTAT=status)
                if ( status .ne. 0 ) then
                   write (UNIT=text, FMT='(I5)') status
                   call mpp_error(FATAL, 'READ_ASCII_FILE: Unable to re-read file '//trim(FILENAME)//'. (IOSTAT = '&
                        //trim(text)//'.')
                else
                   ! A second 'sanity' check on the file
                   pnum_lines = 1

                   do
                      read (UNIT=f_unit, FMT='(A)', IOSTAT=status) str_tmp

                      if ( status .lt. 0 ) then
                         ! deprecate pnum_lines by 1 and ensure pnum_lines is at least 1
                         pnum_lines = max (pnum_lines - 1, 1)
                         exit
                      endif
                      if ( status .gt. 0 ) then
                         write (UNIT=text, FMT='(I5)') pnum_lines
                         call mpp_error(FATAL, 'READ_ASCII_FILE: Error reading line '// &
                                        & trim(text)//' in file '//trim(FILENAME)//'.')
                      end if
                      if(pnum_lines > num_lines) then
                         call mpp_error(FATAL, 'READ_ASCII_FILE: number of lines in file '//trim(FILENAME)// &
                                ' is greater than size(Content(:)). ')
                      end if
                      if ( len_trim(str_tmp) == LENGTH ) then
                         write(UNIT=text, FMT='(I5)') length
                         call mpp_error(FATAL, 'READ_ASCII_FILE: Length of output string ('//trim(text)// &
                                             & ' is too small. Increase the LENGTH value.')
                      end if
                      Content(pnum_lines) = str_tmp
                      pnum_lines = pnum_lines + 1
                   end do
                   if(num_lines .NE. pnum_lines) then
                      call mpp_error(FATAL, 'READ_ASCII_FILE: number of lines in file '//trim(FILENAME)// &
                          ' does not equal to size(Content(:)) ' )
                   end if
                end if
             end if
             close(UNIT=f_unit)
          end if
       else
          call mpp_error(FATAL, 'READ_ASCII_FILE: File '//trim(FILENAME)//' does not exist.')
       end if
    end if

    ! Broadcast character array
    call mpp_broadcast(Content, LENGTH, from_pe, PELIST=PELIST)

  end subroutine read_ascii_file

  !> @brief Produce an inverse permutation. For example, transform [2, 3, 1, 4] to [3, 1, 2, 4].
  !!
  !! The purpose of this subroutine is to convert between (memory dimension) -> (logical dimension) and
  !! (logical dimension) -> (memory dimension) maps.
  !!
  !! @param [in] <x> The original permutation vector
  !! @param [out] <y> The inverted permutation vector
  subroutine inverse_permutation(x, y)
    integer, intent(in) :: x(:)
    integer, intent(out) :: y(size(x))
    integer :: i, n

    y = 0

    n = size(x)
    do i=1,n
      if (x(i).ge.1 .and. x(i).le.n) then
        y(x(i)) = i
      else
        block
          character(2) :: nstr

          write (nstr, "(I0)") n
          call mpp_error(FATAL, "inverse_permutation: Invalid dimension map. &
                                 Values must be in the range from 1 to " // trim(nstr) // ".")
        end block
      endif
    enddo

    if (any(y.eq.0)) then
      call mpp_error(FATAL, "inverse_permutation: Invalid dim_order. Values must be non-repeating.")
    endif
  end subroutine inverse_permutation
# 1398 "mpp/mpp.F90" 2

# 1 "mpp/include/mpp_comm.inc" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

subroutine mpp_init_legacy( flags, localcomm, test_level, alt_input_nml_path )
  integer, optional, intent(in) :: flags !< Flags for debug output, can be MPP_VERBOSE or MPP_DEBUG
  integer, intent(in) :: localcomm !< MPI communicator to use. Only relevant if MPI has already
                                   !! been initialized by an external call to mpi_init.
  integer, optional, intent(in) :: test_level !< Used to exit initialization at certain stages
                                              !! before completion for testing purposes
  character(len=*), optional, intent(in) :: alt_input_nml_path !< Input path for namelist
  type(mpi_comm) :: comm

  comm%mpi_val = localcomm
  call mpp_init_f08(flags, comm, test_level, alt_input_nml_path)
end subroutine mpp_init_legacy





# 1 "mpp/include/mpp_comm_nocomm.inc" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!       ROUTINES TO INITIALIZE/FINALIZE MPP MODULE: mpp_init, mpp_exit        !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> @brief Initialize the @ref mpp_mod module
subroutine mpp_init_f08( flags, localcomm, test_level, alt_input_nml_path )
  integer, optional, intent(in) :: flags !< Flags for debug output, can be MPP_VERBOSE or MPP_DEBUG
  type(mpi_comm), optional, intent(in) :: localcomm !< Id of MPI communicator used to initialize
  integer, optional, intent(in) :: test_level !< Used to exit initialization at certain stages
                                              !! before completion for testing purposes
  character(len=*), optional, intent(in) :: alt_input_nml_path !< Input path for namelist
  integer                       :: my_pe, num_pes, len, i, logunit
  logical                       :: opened, existed
  integer                       :: io_status
  integer                       :: t_level

  if( module_is_initialized )return

  module_is_initialized = .TRUE.
  if(present(test_level)) then
    t_level = test_level
  else
    t_level = -1
  endif
  if(t_level == 0) return

  allocate(peset(0:0))
  !PEsets: make defaults illegal
  peset(:)%count = -1
  peset(:)%comm = MPI_COMM_NULL
  peset(:)%group = MPI_GROUP_NULL
  !0=single-PE, initialized so that count returns 1
  peset(0)%count = 1
  allocate( peset(0)%list(1) )
  peset(0)%list = pe
  current_peset_num = 0
  peset(0)%comm = MPI_COMM_NULL
  world_peset_num = 0
  current_peset_num = world_peset_num !initialize current PEset to world
  if(t_level == 1) return

  !initialize clocks
  call system_clock_default( count=tick0, count_rate=ticks_per_sec, count_max=max_ticks )
  tick_rate = 1./ticks_per_sec
  clock0 = mpp_clock_id( 'Total runtime', flags=MPP_CLOCK_SYNC )
  if(t_level == 2) return

  ! Initialize mpp_datatypes
  ! NOTE: mpp_datatypes is unused in serial mode; this is an empty list
  datatypes%head => null()
  datatypes%tail => null()
  datatypes%length = 0

  ! Create the bytestream (default) mpp_datatype
  ! NOTE: mpp_byte is unused in serial mode
  mpp_byte%counter = -1
  mpp_byte%ndims = -1
  allocate(mpp_byte%sizes(0))
  allocate(mpp_byte%subsizes(0))
  allocate(mpp_byte%starts(0))
  mpp_byte%etype = MPI_DATATYPE_NULL
  mpp_byte%id = MPI_DATATYPE_NULL

  mpp_byte%prev => null()
  mpp_byte%next => null()

  if( PRESENT(flags) )then
     debug   = flags.EQ.MPP_DEBUG
     verbose = flags.EQ.MPP_VERBOSE .OR. debug
  end if
  if(t_level == 3) return

  call mpp_init_logfile()
  if (present(alt_input_nml_path)) then
     call read_input_nml(alt_input_nml_path=alt_input_nml_path)
  else
     call read_input_nml
  end if
  if(t_level == 4) return

  !--- read namelist
  read (input_nml_file, mpp_nml, iostat=io_status)
  if (io_status > 0) then
     call mpp_error(FATAL,'=>mpp_init: Error reading mpp_nml')
  endif
  if(t_level == 5) return

! non-root pe messages written to other location than stdout()
  if (trim(etcfile) /= '/dev/null') then
    write( etcfile,'(a,i6.6)' )trim(etcfile)//'.', pe
  endif
  inquire(file=etcfile, exist=existed)
  if(existed) then
     open( newunit=etc_unit, file=trim(etcfile), status='REPLACE' )
  else
     open( newunit=etc_unit, file=trim(etcfile) )
  endif

  !messages
  if( verbose )call mpp_error( NOTE, 'MPP_INIT: initializing MPP module...' )
  if( pe.EQ.root_pe )then
     logunit = stdlog()
     write( logunit,'(/a)' )'MPP module '//trim(version)
     write( logunit,'(a,i6)' )'MPP started with NPES=', npes
     write( logunit,'(a)' )'Using no library for message passing...'
     write( logunit, '(a26,es12.4,a6,i10,a11)' ) &
          'Realtime clock resolution=', tick_rate, ' sec (', ticks_per_sec, ' ticks/sec)'
     write( logunit, '(a23,es12.4,a6,i20,a7)' ) &
          'Clock rolls over after ', max_ticks*tick_rate, ' sec (', max_ticks, ' ticks)'
  end if

  call mpp_clock_begin(clock0)

  return
end subroutine mpp_init_f08

!#######################################################################
!> @brief To be called at the end of a run
subroutine mpp_exit()
  integer :: i, j, k, n, nmax, istat, out_unit
  real    :: t, tmin, tmax, tavg, tstd
  real    :: m, mmin, mmax, mavg, mstd, t_total
  logical :: opened

  if( .NOT.module_is_initialized )return
  call mpp_set_current_pelist()
  call mpp_clock_end(clock0)
  t_total = clocks(clock0)%total_ticks*tick_rate
  out_unit = stdout()
  if( clock_num.GT.0 )then
     if( ANY(clocks(1:clock_num)%detailed) )then
        call sum_clock_data; call dump_clock_summary
     end if
     if( pe.EQ.root_pe )then
        write( out_unit,'(/a,i6,a)' ) 'Tabulating mpp_clock statistics across ', npes, ' PEs...'
        if( ANY(clocks(1:clock_num)%detailed) ) &
             write( out_unit,'(a)' )'   ... see mpp_clock.out.#### for details on individual PEs.'
        write( out_unit,'(/32x,a)' ) &
             & '      hits          tmin          tmax          tavg          tstd  tfrac grain pemin pemax'
     else
        write( out_unit,'(/37x,a)' ) 'time'
     end if
     call FLUSH( out_unit )
     call mpp_sync()
     do i = 1,clock_num
        if( .NOT.ANY(peset(clocks(i)%peset_num)%list(:).EQ.pe) )cycle
        call mpp_set_current_pelist( peset(clocks(i)%peset_num)%list )
        !times between mpp_clock ticks
        t = clocks(i)%total_ticks*tick_rate
        tmin = t; call mpp_min(tmin)
        tmax = t; call mpp_max(tmax)
        tavg = t; call mpp_sum(tavg); tavg = tavg/mpp_npes()
        tstd = (t-tavg)**2; call mpp_sum(tstd); tstd = sqrt( tstd/mpp_npes() )
        if( pe.EQ.root_pe )write( out_unit,'(a32,i10,4f14.6,f7.3,3i6)' ) &
             clocks(i)%name, clocks(i)%hits, tmin, tmax, tavg, tstd, tavg/t_total, &
             clocks(i)%grain, minval(peset(clocks(i)%peset_num)%list), &
             maxval(peset(clocks(i)%peset_num)%list)
        if (pe.NE.root_pe) write(out_unit,'(a32,f14.6)') clocks(i)%name, clocks(i)%total_ticks*tick_rate
     end do
     if( ANY(clocks(1:clock_num)%detailed) .AND. pe.EQ.root_pe )write( out_unit,'(/32x,a)' ) &
          '       tmin       tmax       tavg       tstd       mmin       mmax       mavg       mstd  mavg/tavg'

     do i = 1,clock_num
        !messages: bytelengths and times
        if( .NOT.clocks(i)%detailed )cycle
        do j = 1,MAX_EVENT_TYPES
           n = clocks(i)%events(j)%calls; nmax = n
           call mpp_max(nmax)
           if( nmax.NE.0 )then
              !don't divide by n because n might be 0
              m = 0
              if( n.GT.0 )m = sum(clocks(i)%events(j)%bytes(1:n))
              mmin = m; call mpp_min(mmin)
              mmax = m; call mpp_max(mmax)
              mavg = m; call mpp_sum(mavg); mavg = mavg/mpp_npes()
              mstd = (m-mavg)**2; call mpp_sum(mstd); mstd = sqrt( mstd/mpp_npes() )
              t = 0
              if( n.GT.0 )t = sum(clocks(i)%events(j)%ticks(1:n))*tick_rate
              tmin = t; call mpp_min(tmin)
              tmax = t; call mpp_max(tmax)
              tavg = t; call mpp_sum(tavg); tavg = tavg/mpp_npes()
              tstd = (t-tavg)**2; call mpp_sum(tstd); tstd = sqrt( tstd/mpp_npes() )
              if( pe.EQ.root_pe )write( out_unit,'(a32,4f11.3,5es11.3)' ) &
                   trim(clocks(i)%name)//' '//trim(clocks(i)%events(j)%name), &
                   tmin, tmax, tavg, tstd, mmin, mmax, mavg, mstd, mavg/tavg
           end if
        end do
     end do
  end if

  inquire(unit=etc_unit, opened=opened)
  if (opened) then
   call FLUSH (etc_unit)
   close(etc_unit)
  endif

  call mpp_set_current_pelist()
  call mpp_sync()
  call mpp_max(mpp_stack_hwm)
  if( pe.EQ.root_pe )write( out_unit,* )'MPP_STACK high water mark=', mpp_stack_hwm

  return
end subroutine mpp_exit

!#######################################################################
  !> Set the mpp_stack variable to be at least n LONG words long
  subroutine mpp_set_stack_size(n)
    integer, intent(in) :: n
    character(len=8)    :: text

    if( n.GT.mpp_stack_size .AND. allocated(mpp_stack) )deallocate(mpp_stack)
    if( .NOT.allocated(mpp_stack) )then
       allocate( mpp_stack(n) )
       mpp_stack_size = n
    end if

    write( text,'(i8)' )n
    if( pe.EQ.root_pe )call mpp_error( NOTE, 'MPP_SET_STACK_SIZE: stack size set to '//text//'.' )

    return
  end subroutine mpp_set_stack_size

    subroutine mpp_broadcast_char(char_data, length, from_pe, pelist )
      character(len=*), intent(inout) :: char_data(:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'mpp_broadcast_text: You must first call mpp_init.' )
      return
    end subroutine mpp_broadcast_char


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                BASIC MESSAGE PASSING ROUTINE: mpp_transmit                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set init value for mpp_type


# 324 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_transmit_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                  MPP_TRANSMIT                               !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> A message-passing routine intended to be reminiscent equally of both MPI and SHMEM
!! put_data and get_data are contiguous real(r8_kind) arrays
!!at each call, your put_data array is put to   to_pe's get_data
!!              your get_data array is got from from_pe's put_data
!!i.e we assume that typically (e.g updating halo regions) each PE performs a put _and_ a get
!!special PE designations:
!!      NULL_PE: to disable a put or a get (e.g at boundaries)
!!      ANY_PE:  if remote PE for the put or get is to be unspecific
!!      ALL_PES: broadcast and collect operations (collect not yet implemented)
!!ideally we would not pass length, but this f77-style call performs better
!!(arrays passed by address, not descriptor) further, this permits <length> contiguous
!!words from an array of any rank to be passed (avoiding f90 rank conformance check)
!!caller is responsible for completion checks (mpp_sync_self) before and after
    subroutine mpp_transmit_real8( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, recv_request, &
                            &  send_request, omp_offload )

      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      real(r8_kind), intent(in)  :: put_data(*)
      real(r8_kind), intent(out) :: get_data(*)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      logical, intent(in), optional :: omp_offload
        ! NOTE: omp_offload is unused in this function

      integer :: i, outunit
      real(r8_kind), allocatable, save :: local_data(:) !local copy used by non-parallel code (no SHMEM or MPI)
      integer(i8_kind),     save :: get_data_addr=-9999
      real(r8_kind)                    :: get_data_local(get_len_nocomm)
      pointer(ptr_get_data, get_data_local)


      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_TRANSMIT: You must first call mpp_init.' )
      if( to_pe.EQ.NULL_PE .AND. from_pe.EQ.NULL_PE )return


      outunit = stdout()
      if( debug )then
          call system_clock_default(tick)
          write( outunit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT begin: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, &
                       &  put_len, get_len
      end if

!do put first and then get
      if( to_pe.GE.0 .AND. to_pe.LT.npes )then
          if( allocated(local_data) ) &
               call mpp_error( FATAL, 'MPP_TRANSMIT: local_data should have been deallocated by prior receive.' )
          if( get_len_nocomm > 0) then  ! pre-post recv
             ptr_get_data = get_data_addr
             do i = 1,get_len_nocomm
                get_data_local(i) = put_data(i)
             end do
             get_len_nocomm = 0
             get_data_addr = -9999
          else
             allocate( local_data(put_len) )
             do i = 1,put_len
                local_data(i) = put_data(i)
             end do
          endif
      else if( to_pe.EQ.ALL_PES )then !this is a broadcast from from_pe
          if( from_pe.LT.0 .OR. from_pe.GE.npes ) &
            call mpp_error( FATAL, 'MPP_TRANSMIT: broadcasting from invalid PE.' )
          if( put_len.GT.get_len ) &
            call mpp_error( FATAL, 'MPP_TRANSMIT: size mismatch between put_data and get_data.' )
          if( pe.EQ.from_pe )then
              if( LOC(get_data).NE.LOC(put_data) )then
!dir$ IVDEP
                  do i = 1,get_len
                     get_data(i) = put_data(i)
                  end do
              end if
          end if
          call mpp_broadcast( get_data, get_len, from_pe )
          return

      else if( to_pe.EQ.ANY_PE )then !we don't have a destination to do puts to, so only do gets

      else if( to_pe.NE.NULL_PE )then  !no other valid cases except NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid to_pe.' )
      end if

!do the get
      if( from_pe.GE.0 .AND. from_pe.LT.npes )then
          if( .NOT.allocated(local_data) ) then
             get_data_addr = LOC(get_data)
             get_len_nocomm  = get_len
          else
             do i = 1,get_len
                get_data(i) = local_data(i)
             end do
             deallocate(local_data)
          endif
      else if( from_pe.EQ.ANY_PE )then

      else if( from_pe.EQ.ALL_PES )then
          call mpp_error( FATAL, 'MPP_TRANSMIT: from_pe=ALL_PES has ambiguous meaning,' &
          & // 'and hence is not implemented.' )

      else if( from_pe.NE.NULL_PE )then !only remaining valid choice is NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid from_pe.' )
      end if

      if( debug )then
          call system_clock_default(tick)
          write( outunit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT end: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, &
                       &  put_len, get_len
      end if
      return
    end subroutine mpp_transmit_real8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_BROADCAST                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_broadcast_real8( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      real(r8_kind), intent(inout) :: data(*)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      integer :: n

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_BROADCAST: You must first call mpp_init.' )
      return
    end subroutine mpp_broadcast_real8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_SCATTER                                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_scatterv_real8( send_data, send_counts, displs, recv_data, recv_count, root_pe, pelist, ierr)
      real(r8_kind), dimension(:),      intent(in) :: send_data
      real(r8_kind), dimension(:,:,:),  intent(inout) :: recv_data
      integer,                      intent(in) :: recv_count, root_pe
      integer, dimension(:),        intent(in) :: send_counts, displs, pelist
      integer,                      intent(inout) :: ierr

      integer :: n, i, j, k

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SCATTERV: You must first call mpp_init.' )

      n = 1
      do k = 1, size(recv_data, 3)
        do j = 1, size(recv_data, 2)
          do i = 1, size(recv_data, 1)
            recv_data(i,j,k) = send_data(n)
            n = n + 1
          enddo
        enddo
      enddo

    end subroutine mpp_scatterv_real8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_GATHER                                   !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_gather_real8( send_data, recv_data, count, root_pe, pelist, ierr)
      real(r8_kind), dimension(:), intent(in) :: send_data
      real(r8_kind), dimension(:), intent(inout) :: recv_data
      integer, dimension(:), intent(in) :: pelist(:)
      integer, intent(in) :: count, root_pe
      integer, intent(inout) :: ierr

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GATHER: You must first call mpp_init.' )

      recv_data = send_data

    end subroutine mpp_gather_real8

    subroutine mpp_gatherv_real8( send_data, send_count, recv_data, recv_counts, displs, root_pe, pelist, ierr)
      real(r8_kind), dimension(:), intent(in) :: send_data
      real(r8_kind), dimension(:), intent(inout) :: recv_data
      integer, dimension(:), intent(in) :: recv_counts, displs, pelist
      integer, intent(in) :: send_count, root_pe
      integer, intent(inout) :: ierr

      if( .NOT.module_is_initialized ) call mpp_error( FATAL, 'MPP_GATHERV: You must first call mpp_init.' )

      recv_data = send_data

    end subroutine mpp_gatherv_real8

!####################################################################################

# 1 "mpp/include/mpp_transmit.inc" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                  MPP_TRANSMIT                               !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_transmit_real8_scalar( put_data, to_pe, get_data, from_pe, plen, glen, block, tag, &
                                    recv_request, send_request)
      integer, intent(in) :: to_pe, from_pe
      real(r8_kind), intent(in)  :: put_data
      real(r8_kind), intent(out) :: get_data
      integer, optional,  intent(in) :: plen, glen
      logical, intent(in),  optional :: block
      integer, intent(in),  optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      integer                       :: put_len, get_len
      real(r8_kind) :: put_data1D(1), get_data1D(1)
      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )

      get_data = 0.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      put_len=1; if(PRESENT(plen))put_len=plen
      get_len=1; if(PRESENT(glen))get_len=glen
      call mpp_transmit_real8 ( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                           recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_real8_scalar

    subroutine mpp_transmit_real8_2d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      real(r8_kind), intent(in)  :: put_data(:,:)
      real(r8_kind), intent(out) :: get_data(:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      real(r8_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = 0.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_real8_2d

    subroutine mpp_transmit_real8_3d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      real(r8_kind), intent(in)  :: put_data(:,:,:)
      real(r8_kind), intent(out) :: get_data(:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      real(r8_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = 0.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_real8_3d

    subroutine mpp_transmit_real8_4d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      real(r8_kind), intent(in)  :: put_data(:,:,:,:)
      real(r8_kind), intent(out) :: get_data(:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      real(r8_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = 0.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_real8_4d

    subroutine mpp_transmit_real8_5d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      real(r8_kind), intent(in)  :: put_data(:,:,:,:,:)
      real(r8_kind), intent(out) :: get_data(:,:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      real(r8_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = 0.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_real8_5d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                              MPP_SEND and RECV                              !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_recv_real8( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      real(r8_kind), intent(out) :: get_data(*)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      real(r8_kind) :: dummy(1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_real8

    subroutine mpp_send_real8( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      real(r8_kind), intent(in) :: put_data(*)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      real(r8_kind) :: dummy(1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag=tag, send_request=request )
    end subroutine mpp_send_real8

    subroutine mpp_recv_real8_scalar( get_data, from_pe, glen, block, tag, request, omp_offload )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: from_pe
      real(r8_kind), intent(out) :: get_data
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      integer, optional, intent(in) :: glen
      logical, optional, intent(in) :: omp_offload
      integer                       :: get_len
      real(r8_kind) :: get_data1D(1)
      real(r8_kind) :: dummy(1)

      pointer( ptr, get_data1D )
      get_data = 0.

      ptr = LOC(get_data)
      get_len=1; if(PRESENT(glen))get_len=glen
      call mpp_transmit( dummy, 1, NULL_PE, get_data1D, get_len, from_pe, &
          block, tag, recv_request=request, omp_offload=omp_offload )

    end subroutine mpp_recv_real8_scalar

    subroutine mpp_send_real8_scalar( put_data, to_pe, plen, tag, request, omp_offload)
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: to_pe
      real(r8_kind), intent(in) :: put_data
      integer, optional, intent(in) :: plen
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      logical, optional, intent(in) :: omp_offload
      integer                       :: put_len
      real(r8_kind) :: put_data1D(1)
      real(r8_kind) :: dummy(1)

      pointer( ptr, put_data1D )
      ptr = LOC(put_data)
      put_len=1; if(PRESENT(plen))put_len=plen
      call mpp_transmit( put_data1D, put_len, to_pe, dummy, 1, NULL_PE, &
          tag=tag, send_request=request, omp_offload=omp_offload )

    end subroutine mpp_send_real8_scalar

    subroutine mpp_recv_real8_2d( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      real(r8_kind), intent(out) :: get_data(:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      real(r8_kind) :: dummy(1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, &
          block, tag, recv_request=request )
    end subroutine mpp_recv_real8_2d

    subroutine mpp_send_real8_2d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      real(r8_kind), intent(in) :: put_data(:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      real(r8_kind) :: dummy(1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_real8_2d

    subroutine mpp_recv_real8_3d( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      real(r8_kind), intent(out) :: get_data(:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      real(r8_kind) :: dummy(1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_real8_3d

    subroutine mpp_send_real8_3d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      real(r8_kind), intent(in) :: put_data(:,:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      real(r8_kind) :: dummy(1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_real8_3d

    subroutine mpp_recv_real8_4d( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      real(r8_kind), intent(out) :: get_data(:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      real(r8_kind) :: dummy(1,1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_real8_4d

    subroutine mpp_send_real8_4d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      real(r8_kind), intent(in) :: put_data(:,:,:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      real(r8_kind) :: dummy(1,1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_real8_4d

    subroutine mpp_recv_real8_5d( get_data, get_len, from_pe, block, tag, request)
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      real(r8_kind), intent(out) :: get_data(:,:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      real(r8_kind) :: dummy(1,1,1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_real8_5d

    subroutine mpp_send_real8_5d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      real(r8_kind), intent(in) :: put_data(:,:,:,:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      real(r8_kind) :: dummy(1,1,1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_real8_5d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_BROADCAST                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_broadcast_real8_scalar( broadcast_data, from_pe, pelist )
      real(r8_kind), intent(inout) :: broadcast_data
      integer, intent(in) :: from_pe
      integer, intent(in), optional :: pelist(:)
      real(r8_kind) :: data1D(1)

      pointer( ptr, data1D )

      ptr = LOC(broadcast_data)
      call mpp_broadcast_real8( data1D, 1, from_pe, pelist )

      return
    end subroutine mpp_broadcast_real8_scalar

    subroutine mpp_broadcast_real8_2d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      real(r8_kind), intent(inout) :: broadcast_data(:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      real(r8_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine mpp_broadcast_real8_2d

    subroutine mpp_broadcast_real8_3d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      real(r8_kind), intent(inout) :: broadcast_data(:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      real(r8_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
   end subroutine mpp_broadcast_real8_3d

    subroutine mpp_broadcast_real8_4d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      real(r8_kind), intent(inout) :: broadcast_data(:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      real(r8_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine mpp_broadcast_real8_4d

    subroutine mpp_broadcast_real8_5d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      real(r8_kind), intent(inout) :: broadcast_data(:,:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      real(r8_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine mpp_broadcast_real8_5d
# 217 "mpp/include/mpp_transmit_nocomm.fh" 2
# 324 "mpp/include/mpp_comm_nocomm.inc" 2

# 388 "mpp/include/mpp_comm_nocomm.inc"

# 450 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_transmit_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                  MPP_TRANSMIT                               !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> A message-passing routine intended to be reminiscent equally of both MPI and SHMEM
!! put_data and get_data are contiguous real(r4_kind) arrays
!!at each call, your put_data array is put to   to_pe's get_data
!!              your get_data array is got from from_pe's put_data
!!i.e we assume that typically (e.g updating halo regions) each PE performs a put _and_ a get
!!special PE designations:
!!      NULL_PE: to disable a put or a get (e.g at boundaries)
!!      ANY_PE:  if remote PE for the put or get is to be unspecific
!!      ALL_PES: broadcast and collect operations (collect not yet implemented)
!!ideally we would not pass length, but this f77-style call performs better
!!(arrays passed by address, not descriptor) further, this permits <length> contiguous
!!words from an array of any rank to be passed (avoiding f90 rank conformance check)
!!caller is responsible for completion checks (mpp_sync_self) before and after
    subroutine mpp_transmit_real4( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, recv_request, &
                            &  send_request, omp_offload )

      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      real(r4_kind), intent(in)  :: put_data(*)
      real(r4_kind), intent(out) :: get_data(*)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      logical, intent(in), optional :: omp_offload
        ! NOTE: omp_offload is unused in this function

      integer :: i, outunit
      real(r4_kind), allocatable, save :: local_data(:) !local copy used by non-parallel code (no SHMEM or MPI)
      integer(i8_kind),     save :: get_data_addr=-9999
      real(r4_kind)                    :: get_data_local(get_len_nocomm)
      pointer(ptr_get_data, get_data_local)


      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_TRANSMIT: You must first call mpp_init.' )
      if( to_pe.EQ.NULL_PE .AND. from_pe.EQ.NULL_PE )return


      outunit = stdout()
      if( debug )then
          call system_clock_default(tick)
          write( outunit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT begin: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, &
                       &  put_len, get_len
      end if

!do put first and then get
      if( to_pe.GE.0 .AND. to_pe.LT.npes )then
          if( allocated(local_data) ) &
               call mpp_error( FATAL, 'MPP_TRANSMIT: local_data should have been deallocated by prior receive.' )
          if( get_len_nocomm > 0) then  ! pre-post recv
             ptr_get_data = get_data_addr
             do i = 1,get_len_nocomm
                get_data_local(i) = put_data(i)
             end do
             get_len_nocomm = 0
             get_data_addr = -9999
          else
             allocate( local_data(put_len) )
             do i = 1,put_len
                local_data(i) = put_data(i)
             end do
          endif
      else if( to_pe.EQ.ALL_PES )then !this is a broadcast from from_pe
          if( from_pe.LT.0 .OR. from_pe.GE.npes ) &
            call mpp_error( FATAL, 'MPP_TRANSMIT: broadcasting from invalid PE.' )
          if( put_len.GT.get_len ) &
            call mpp_error( FATAL, 'MPP_TRANSMIT: size mismatch between put_data and get_data.' )
          if( pe.EQ.from_pe )then
              if( LOC(get_data).NE.LOC(put_data) )then
!dir$ IVDEP
                  do i = 1,get_len
                     get_data(i) = put_data(i)
                  end do
              end if
          end if
          call mpp_broadcast( get_data, get_len, from_pe )
          return

      else if( to_pe.EQ.ANY_PE )then !we don't have a destination to do puts to, so only do gets

      else if( to_pe.NE.NULL_PE )then  !no other valid cases except NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid to_pe.' )
      end if

!do the get
      if( from_pe.GE.0 .AND. from_pe.LT.npes )then
          if( .NOT.allocated(local_data) ) then
             get_data_addr = LOC(get_data)
             get_len_nocomm  = get_len
          else
             do i = 1,get_len
                get_data(i) = local_data(i)
             end do
             deallocate(local_data)
          endif
      else if( from_pe.EQ.ANY_PE )then

      else if( from_pe.EQ.ALL_PES )then
          call mpp_error( FATAL, 'MPP_TRANSMIT: from_pe=ALL_PES has ambiguous meaning,' &
          & // 'and hence is not implemented.' )

      else if( from_pe.NE.NULL_PE )then !only remaining valid choice is NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid from_pe.' )
      end if

      if( debug )then
          call system_clock_default(tick)
          write( outunit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT end: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, &
                       &  put_len, get_len
      end if
      return
    end subroutine mpp_transmit_real4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_BROADCAST                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_broadcast_real4( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      real(r4_kind), intent(inout) :: data(*)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      integer :: n

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_BROADCAST: You must first call mpp_init.' )
      return
    end subroutine mpp_broadcast_real4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_SCATTER                                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_scatterv_real4( send_data, send_counts, displs, recv_data, recv_count, root_pe, pelist, ierr)
      real(r4_kind), dimension(:),      intent(in) :: send_data
      real(r4_kind), dimension(:,:,:),  intent(inout) :: recv_data
      integer,                      intent(in) :: recv_count, root_pe
      integer, dimension(:),        intent(in) :: send_counts, displs, pelist
      integer,                      intent(inout) :: ierr

      integer :: n, i, j, k

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SCATTERV: You must first call mpp_init.' )

      n = 1
      do k = 1, size(recv_data, 3)
        do j = 1, size(recv_data, 2)
          do i = 1, size(recv_data, 1)
            recv_data(i,j,k) = send_data(n)
            n = n + 1
          enddo
        enddo
      enddo

    end subroutine mpp_scatterv_real4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_GATHER                                   !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_gather_real4( send_data, recv_data, count, root_pe, pelist, ierr)
      real(r4_kind), dimension(:), intent(in) :: send_data
      real(r4_kind), dimension(:), intent(inout) :: recv_data
      integer, dimension(:), intent(in) :: pelist(:)
      integer, intent(in) :: count, root_pe
      integer, intent(inout) :: ierr

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GATHER: You must first call mpp_init.' )

      recv_data = send_data

    end subroutine mpp_gather_real4

    subroutine mpp_gatherv_real4( send_data, send_count, recv_data, recv_counts, displs, root_pe, pelist, ierr)
      real(r4_kind), dimension(:), intent(in) :: send_data
      real(r4_kind), dimension(:), intent(inout) :: recv_data
      integer, dimension(:), intent(in) :: recv_counts, displs, pelist
      integer, intent(in) :: send_count, root_pe
      integer, intent(inout) :: ierr

      if( .NOT.module_is_initialized ) call mpp_error( FATAL, 'MPP_GATHERV: You must first call mpp_init.' )

      recv_data = send_data

    end subroutine mpp_gatherv_real4

!####################################################################################

# 1 "mpp/include/mpp_transmit.inc" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                  MPP_TRANSMIT                               !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_transmit_real4_scalar( put_data, to_pe, get_data, from_pe, plen, glen, block, tag, &
                                    recv_request, send_request)
      integer, intent(in) :: to_pe, from_pe
      real(r4_kind), intent(in)  :: put_data
      real(r4_kind), intent(out) :: get_data
      integer, optional,  intent(in) :: plen, glen
      logical, intent(in),  optional :: block
      integer, intent(in),  optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      integer                       :: put_len, get_len
      real(r4_kind) :: put_data1D(1), get_data1D(1)
      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )

      get_data = 0.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      put_len=1; if(PRESENT(plen))put_len=plen
      get_len=1; if(PRESENT(glen))get_len=glen
      call mpp_transmit_real4 ( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                           recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_real4_scalar

    subroutine mpp_transmit_real4_2d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      real(r4_kind), intent(in)  :: put_data(:,:)
      real(r4_kind), intent(out) :: get_data(:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      real(r4_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = 0.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_real4_2d

    subroutine mpp_transmit_real4_3d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      real(r4_kind), intent(in)  :: put_data(:,:,:)
      real(r4_kind), intent(out) :: get_data(:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      real(r4_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = 0.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_real4_3d

    subroutine mpp_transmit_real4_4d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      real(r4_kind), intent(in)  :: put_data(:,:,:,:)
      real(r4_kind), intent(out) :: get_data(:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      real(r4_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = 0.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_real4_4d

    subroutine mpp_transmit_real4_5d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      real(r4_kind), intent(in)  :: put_data(:,:,:,:,:)
      real(r4_kind), intent(out) :: get_data(:,:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      real(r4_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = 0.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_real4_5d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                              MPP_SEND and RECV                              !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_recv_real4( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      real(r4_kind), intent(out) :: get_data(*)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      real(r4_kind) :: dummy(1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_real4

    subroutine mpp_send_real4( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      real(r4_kind), intent(in) :: put_data(*)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      real(r4_kind) :: dummy(1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag=tag, send_request=request )
    end subroutine mpp_send_real4

    subroutine mpp_recv_real4_scalar( get_data, from_pe, glen, block, tag, request, omp_offload )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: from_pe
      real(r4_kind), intent(out) :: get_data
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      integer, optional, intent(in) :: glen
      logical, optional, intent(in) :: omp_offload
      integer                       :: get_len
      real(r4_kind) :: get_data1D(1)
      real(r4_kind) :: dummy(1)

      pointer( ptr, get_data1D )
      get_data = 0.

      ptr = LOC(get_data)
      get_len=1; if(PRESENT(glen))get_len=glen
      call mpp_transmit( dummy, 1, NULL_PE, get_data1D, get_len, from_pe, &
          block, tag, recv_request=request, omp_offload=omp_offload )

    end subroutine mpp_recv_real4_scalar

    subroutine mpp_send_real4_scalar( put_data, to_pe, plen, tag, request, omp_offload)
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: to_pe
      real(r4_kind), intent(in) :: put_data
      integer, optional, intent(in) :: plen
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      logical, optional, intent(in) :: omp_offload
      integer                       :: put_len
      real(r4_kind) :: put_data1D(1)
      real(r4_kind) :: dummy(1)

      pointer( ptr, put_data1D )
      ptr = LOC(put_data)
      put_len=1; if(PRESENT(plen))put_len=plen
      call mpp_transmit( put_data1D, put_len, to_pe, dummy, 1, NULL_PE, &
          tag=tag, send_request=request, omp_offload=omp_offload )

    end subroutine mpp_send_real4_scalar

    subroutine mpp_recv_real4_2d( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      real(r4_kind), intent(out) :: get_data(:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      real(r4_kind) :: dummy(1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, &
          block, tag, recv_request=request )
    end subroutine mpp_recv_real4_2d

    subroutine mpp_send_real4_2d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      real(r4_kind), intent(in) :: put_data(:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      real(r4_kind) :: dummy(1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_real4_2d

    subroutine mpp_recv_real4_3d( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      real(r4_kind), intent(out) :: get_data(:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      real(r4_kind) :: dummy(1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_real4_3d

    subroutine mpp_send_real4_3d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      real(r4_kind), intent(in) :: put_data(:,:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      real(r4_kind) :: dummy(1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_real4_3d

    subroutine mpp_recv_real4_4d( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      real(r4_kind), intent(out) :: get_data(:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      real(r4_kind) :: dummy(1,1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_real4_4d

    subroutine mpp_send_real4_4d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      real(r4_kind), intent(in) :: put_data(:,:,:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      real(r4_kind) :: dummy(1,1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_real4_4d

    subroutine mpp_recv_real4_5d( get_data, get_len, from_pe, block, tag, request)
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      real(r4_kind), intent(out) :: get_data(:,:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      real(r4_kind) :: dummy(1,1,1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_real4_5d

    subroutine mpp_send_real4_5d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      real(r4_kind), intent(in) :: put_data(:,:,:,:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      real(r4_kind) :: dummy(1,1,1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_real4_5d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_BROADCAST                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_broadcast_real4_scalar( broadcast_data, from_pe, pelist )
      real(r4_kind), intent(inout) :: broadcast_data
      integer, intent(in) :: from_pe
      integer, intent(in), optional :: pelist(:)
      real(r4_kind) :: data1D(1)

      pointer( ptr, data1D )

      ptr = LOC(broadcast_data)
      call mpp_broadcast_real4( data1D, 1, from_pe, pelist )

      return
    end subroutine mpp_broadcast_real4_scalar

    subroutine mpp_broadcast_real4_2d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      real(r4_kind), intent(inout) :: broadcast_data(:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      real(r4_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine mpp_broadcast_real4_2d

    subroutine mpp_broadcast_real4_3d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      real(r4_kind), intent(inout) :: broadcast_data(:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      real(r4_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
   end subroutine mpp_broadcast_real4_3d

    subroutine mpp_broadcast_real4_4d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      real(r4_kind), intent(inout) :: broadcast_data(:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      real(r4_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine mpp_broadcast_real4_4d

    subroutine mpp_broadcast_real4_5d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      real(r4_kind), intent(inout) :: broadcast_data(:,:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      real(r4_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine mpp_broadcast_real4_5d
# 217 "mpp/include/mpp_transmit_nocomm.fh" 2
# 450 "mpp/include/mpp_comm_nocomm.inc" 2

# 514 "mpp/include/mpp_comm_nocomm.inc"

# 578 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_transmit_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                  MPP_TRANSMIT                               !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> A message-passing routine intended to be reminiscent equally of both MPI and SHMEM
!! put_data and get_data are contiguous integer(i8_kind) arrays
!!at each call, your put_data array is put to   to_pe's get_data
!!              your get_data array is got from from_pe's put_data
!!i.e we assume that typically (e.g updating halo regions) each PE performs a put _and_ a get
!!special PE designations:
!!      NULL_PE: to disable a put or a get (e.g at boundaries)
!!      ANY_PE:  if remote PE for the put or get is to be unspecific
!!      ALL_PES: broadcast and collect operations (collect not yet implemented)
!!ideally we would not pass length, but this f77-style call performs better
!!(arrays passed by address, not descriptor) further, this permits <length> contiguous
!!words from an array of any rank to be passed (avoiding f90 rank conformance check)
!!caller is responsible for completion checks (mpp_sync_self) before and after
    subroutine mpp_transmit_int8( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, recv_request, &
                            &  send_request, omp_offload )

      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      integer(i8_kind), intent(in)  :: put_data(*)
      integer(i8_kind), intent(out) :: get_data(*)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      logical, intent(in), optional :: omp_offload
        ! NOTE: omp_offload is unused in this function

      integer :: i, outunit
      integer(i8_kind), allocatable, save :: local_data(:) !local copy used by non-parallel code (no SHMEM or MPI)
      integer(i8_kind),     save :: get_data_addr=-9999
      integer(i8_kind)                    :: get_data_local(get_len_nocomm)
      pointer(ptr_get_data, get_data_local)


      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_TRANSMIT: You must first call mpp_init.' )
      if( to_pe.EQ.NULL_PE .AND. from_pe.EQ.NULL_PE )return


      outunit = stdout()
      if( debug )then
          call system_clock_default(tick)
          write( outunit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT begin: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, &
                       &  put_len, get_len
      end if

!do put first and then get
      if( to_pe.GE.0 .AND. to_pe.LT.npes )then
          if( allocated(local_data) ) &
               call mpp_error( FATAL, 'MPP_TRANSMIT: local_data should have been deallocated by prior receive.' )
          if( get_len_nocomm > 0) then  ! pre-post recv
             ptr_get_data = get_data_addr
             do i = 1,get_len_nocomm
                get_data_local(i) = put_data(i)
             end do
             get_len_nocomm = 0
             get_data_addr = -9999
          else
             allocate( local_data(put_len) )
             do i = 1,put_len
                local_data(i) = put_data(i)
             end do
          endif
      else if( to_pe.EQ.ALL_PES )then !this is a broadcast from from_pe
          if( from_pe.LT.0 .OR. from_pe.GE.npes ) &
            call mpp_error( FATAL, 'MPP_TRANSMIT: broadcasting from invalid PE.' )
          if( put_len.GT.get_len ) &
            call mpp_error( FATAL, 'MPP_TRANSMIT: size mismatch between put_data and get_data.' )
          if( pe.EQ.from_pe )then
              if( LOC(get_data).NE.LOC(put_data) )then
!dir$ IVDEP
                  do i = 1,get_len
                     get_data(i) = put_data(i)
                  end do
              end if
          end if
          call mpp_broadcast( get_data, get_len, from_pe )
          return

      else if( to_pe.EQ.ANY_PE )then !we don't have a destination to do puts to, so only do gets

      else if( to_pe.NE.NULL_PE )then  !no other valid cases except NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid to_pe.' )
      end if

!do the get
      if( from_pe.GE.0 .AND. from_pe.LT.npes )then
          if( .NOT.allocated(local_data) ) then
             get_data_addr = LOC(get_data)
             get_len_nocomm  = get_len
          else
             do i = 1,get_len
                get_data(i) = local_data(i)
             end do
             deallocate(local_data)
          endif
      else if( from_pe.EQ.ANY_PE )then

      else if( from_pe.EQ.ALL_PES )then
          call mpp_error( FATAL, 'MPP_TRANSMIT: from_pe=ALL_PES has ambiguous meaning,' &
          & // 'and hence is not implemented.' )

      else if( from_pe.NE.NULL_PE )then !only remaining valid choice is NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid from_pe.' )
      end if

      if( debug )then
          call system_clock_default(tick)
          write( outunit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT end: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, &
                       &  put_len, get_len
      end if
      return
    end subroutine mpp_transmit_int8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_BROADCAST                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_broadcast_int8( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      integer(i8_kind), intent(inout) :: data(*)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      integer :: n

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_BROADCAST: You must first call mpp_init.' )
      return
    end subroutine mpp_broadcast_int8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_SCATTER                                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_scatterv_int8( send_data, send_counts, displs, recv_data, recv_count, root_pe, pelist, ierr)
      integer(i8_kind), dimension(:),      intent(in) :: send_data
      integer(i8_kind), dimension(:,:,:),  intent(inout) :: recv_data
      integer,                      intent(in) :: recv_count, root_pe
      integer, dimension(:),        intent(in) :: send_counts, displs, pelist
      integer,                      intent(inout) :: ierr

      integer :: n, i, j, k

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SCATTERV: You must first call mpp_init.' )

      n = 1
      do k = 1, size(recv_data, 3)
        do j = 1, size(recv_data, 2)
          do i = 1, size(recv_data, 1)
            recv_data(i,j,k) = send_data(n)
            n = n + 1
          enddo
        enddo
      enddo

    end subroutine mpp_scatterv_int8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_GATHER                                   !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_gather_int8( send_data, recv_data, count, root_pe, pelist, ierr)
      integer(i8_kind), dimension(:), intent(in) :: send_data
      integer(i8_kind), dimension(:), intent(inout) :: recv_data
      integer, dimension(:), intent(in) :: pelist(:)
      integer, intent(in) :: count, root_pe
      integer, intent(inout) :: ierr

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GATHER: You must first call mpp_init.' )

      recv_data = send_data

    end subroutine mpp_gather_int8

    subroutine mpp_gatherv_int8( send_data, send_count, recv_data, recv_counts, displs, root_pe, pelist, ierr)
      integer(i8_kind), dimension(:), intent(in) :: send_data
      integer(i8_kind), dimension(:), intent(inout) :: recv_data
      integer, dimension(:), intent(in) :: recv_counts, displs, pelist
      integer, intent(in) :: send_count, root_pe
      integer, intent(inout) :: ierr

      if( .NOT.module_is_initialized ) call mpp_error( FATAL, 'MPP_GATHERV: You must first call mpp_init.' )

      recv_data = send_data

    end subroutine mpp_gatherv_int8

!####################################################################################

# 1 "mpp/include/mpp_transmit.inc" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                  MPP_TRANSMIT                               !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_transmit_int8_scalar( put_data, to_pe, get_data, from_pe, plen, glen, block, tag, &
                                    recv_request, send_request)
      integer, intent(in) :: to_pe, from_pe
      integer(i8_kind), intent(in)  :: put_data
      integer(i8_kind), intent(out) :: get_data
      integer, optional,  intent(in) :: plen, glen
      logical, intent(in),  optional :: block
      integer, intent(in),  optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      integer                       :: put_len, get_len
      integer(i8_kind) :: put_data1D(1), get_data1D(1)
      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )

      get_data = 0

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      put_len=1; if(PRESENT(plen))put_len=plen
      get_len=1; if(PRESENT(glen))get_len=glen
      call mpp_transmit_int8 ( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                           recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_int8_scalar

    subroutine mpp_transmit_int8_2d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      integer(i8_kind), intent(in)  :: put_data(:,:)
      integer(i8_kind), intent(out) :: get_data(:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      integer(i8_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = 0

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_int8_2d

    subroutine mpp_transmit_int8_3d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      integer(i8_kind), intent(in)  :: put_data(:,:,:)
      integer(i8_kind), intent(out) :: get_data(:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      integer(i8_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = 0

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_int8_3d

    subroutine mpp_transmit_int8_4d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      integer(i8_kind), intent(in)  :: put_data(:,:,:,:)
      integer(i8_kind), intent(out) :: get_data(:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      integer(i8_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = 0

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_int8_4d

    subroutine mpp_transmit_int8_5d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      integer(i8_kind), intent(in)  :: put_data(:,:,:,:,:)
      integer(i8_kind), intent(out) :: get_data(:,:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      integer(i8_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = 0

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_int8_5d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                              MPP_SEND and RECV                              !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_recv_int8( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      integer(i8_kind), intent(out) :: get_data(*)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      integer(i8_kind) :: dummy(1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_int8

    subroutine mpp_send_int8( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      integer(i8_kind), intent(in) :: put_data(*)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      integer(i8_kind) :: dummy(1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag=tag, send_request=request )
    end subroutine mpp_send_int8

    subroutine mpp_recv_int8_scalar( get_data, from_pe, glen, block, tag, request, omp_offload )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: from_pe
      integer(i8_kind), intent(out) :: get_data
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      integer, optional, intent(in) :: glen
      logical, optional, intent(in) :: omp_offload
      integer                       :: get_len
      integer(i8_kind) :: get_data1D(1)
      integer(i8_kind) :: dummy(1)

      pointer( ptr, get_data1D )
      get_data = 0

      ptr = LOC(get_data)
      get_len=1; if(PRESENT(glen))get_len=glen
      call mpp_transmit( dummy, 1, NULL_PE, get_data1D, get_len, from_pe, &
          block, tag, recv_request=request, omp_offload=omp_offload )

    end subroutine mpp_recv_int8_scalar

    subroutine mpp_send_int8_scalar( put_data, to_pe, plen, tag, request, omp_offload)
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: to_pe
      integer(i8_kind), intent(in) :: put_data
      integer, optional, intent(in) :: plen
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      logical, optional, intent(in) :: omp_offload
      integer                       :: put_len
      integer(i8_kind) :: put_data1D(1)
      integer(i8_kind) :: dummy(1)

      pointer( ptr, put_data1D )
      ptr = LOC(put_data)
      put_len=1; if(PRESENT(plen))put_len=plen
      call mpp_transmit( put_data1D, put_len, to_pe, dummy, 1, NULL_PE, &
          tag=tag, send_request=request, omp_offload=omp_offload )

    end subroutine mpp_send_int8_scalar

    subroutine mpp_recv_int8_2d( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      integer(i8_kind), intent(out) :: get_data(:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      integer(i8_kind) :: dummy(1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, &
          block, tag, recv_request=request )
    end subroutine mpp_recv_int8_2d

    subroutine mpp_send_int8_2d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      integer(i8_kind), intent(in) :: put_data(:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      integer(i8_kind) :: dummy(1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_int8_2d

    subroutine mpp_recv_int8_3d( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      integer(i8_kind), intent(out) :: get_data(:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      integer(i8_kind) :: dummy(1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_int8_3d

    subroutine mpp_send_int8_3d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      integer(i8_kind), intent(in) :: put_data(:,:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      integer(i8_kind) :: dummy(1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_int8_3d

    subroutine mpp_recv_int8_4d( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      integer(i8_kind), intent(out) :: get_data(:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      integer(i8_kind) :: dummy(1,1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_int8_4d

    subroutine mpp_send_int8_4d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      integer(i8_kind), intent(in) :: put_data(:,:,:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      integer(i8_kind) :: dummy(1,1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_int8_4d

    subroutine mpp_recv_int8_5d( get_data, get_len, from_pe, block, tag, request)
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      integer(i8_kind), intent(out) :: get_data(:,:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      integer(i8_kind) :: dummy(1,1,1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_int8_5d

    subroutine mpp_send_int8_5d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      integer(i8_kind), intent(in) :: put_data(:,:,:,:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      integer(i8_kind) :: dummy(1,1,1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_int8_5d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_BROADCAST                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_broadcast_int8_scalar( broadcast_data, from_pe, pelist )
      integer(i8_kind), intent(inout) :: broadcast_data
      integer, intent(in) :: from_pe
      integer, intent(in), optional :: pelist(:)
      integer(i8_kind) :: data1D(1)

      pointer( ptr, data1D )

      ptr = LOC(broadcast_data)
      call mpp_broadcast_int8( data1D, 1, from_pe, pelist )

      return
    end subroutine mpp_broadcast_int8_scalar

    subroutine mpp_broadcast_int8_2d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      integer(i8_kind), intent(inout) :: broadcast_data(:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      integer(i8_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine mpp_broadcast_int8_2d

    subroutine mpp_broadcast_int8_3d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      integer(i8_kind), intent(inout) :: broadcast_data(:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      integer(i8_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
   end subroutine mpp_broadcast_int8_3d

    subroutine mpp_broadcast_int8_4d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      integer(i8_kind), intent(inout) :: broadcast_data(:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      integer(i8_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine mpp_broadcast_int8_4d

    subroutine mpp_broadcast_int8_5d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      integer(i8_kind), intent(inout) :: broadcast_data(:,:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      integer(i8_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine mpp_broadcast_int8_5d
# 217 "mpp/include/mpp_transmit_nocomm.fh" 2
# 578 "mpp/include/mpp_comm_nocomm.inc" 2

# 640 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_transmit_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                  MPP_TRANSMIT                               !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> A message-passing routine intended to be reminiscent equally of both MPI and SHMEM
!! put_data and get_data are contiguous integer(i4_kind) arrays
!!at each call, your put_data array is put to   to_pe's get_data
!!              your get_data array is got from from_pe's put_data
!!i.e we assume that typically (e.g updating halo regions) each PE performs a put _and_ a get
!!special PE designations:
!!      NULL_PE: to disable a put or a get (e.g at boundaries)
!!      ANY_PE:  if remote PE for the put or get is to be unspecific
!!      ALL_PES: broadcast and collect operations (collect not yet implemented)
!!ideally we would not pass length, but this f77-style call performs better
!!(arrays passed by address, not descriptor) further, this permits <length> contiguous
!!words from an array of any rank to be passed (avoiding f90 rank conformance check)
!!caller is responsible for completion checks (mpp_sync_self) before and after
    subroutine mpp_transmit_int4( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, recv_request, &
                            &  send_request, omp_offload )

      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      integer(i4_kind), intent(in)  :: put_data(*)
      integer(i4_kind), intent(out) :: get_data(*)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      logical, intent(in), optional :: omp_offload
        ! NOTE: omp_offload is unused in this function

      integer :: i, outunit
      integer(i4_kind), allocatable, save :: local_data(:) !local copy used by non-parallel code (no SHMEM or MPI)
      integer(i8_kind),     save :: get_data_addr=-9999
      integer(i4_kind)                    :: get_data_local(get_len_nocomm)
      pointer(ptr_get_data, get_data_local)


      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_TRANSMIT: You must first call mpp_init.' )
      if( to_pe.EQ.NULL_PE .AND. from_pe.EQ.NULL_PE )return


      outunit = stdout()
      if( debug )then
          call system_clock_default(tick)
          write( outunit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT begin: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, &
                       &  put_len, get_len
      end if

!do put first and then get
      if( to_pe.GE.0 .AND. to_pe.LT.npes )then
          if( allocated(local_data) ) &
               call mpp_error( FATAL, 'MPP_TRANSMIT: local_data should have been deallocated by prior receive.' )
          if( get_len_nocomm > 0) then  ! pre-post recv
             ptr_get_data = get_data_addr
             do i = 1,get_len_nocomm
                get_data_local(i) = put_data(i)
             end do
             get_len_nocomm = 0
             get_data_addr = -9999
          else
             allocate( local_data(put_len) )
             do i = 1,put_len
                local_data(i) = put_data(i)
             end do
          endif
      else if( to_pe.EQ.ALL_PES )then !this is a broadcast from from_pe
          if( from_pe.LT.0 .OR. from_pe.GE.npes ) &
            call mpp_error( FATAL, 'MPP_TRANSMIT: broadcasting from invalid PE.' )
          if( put_len.GT.get_len ) &
            call mpp_error( FATAL, 'MPP_TRANSMIT: size mismatch between put_data and get_data.' )
          if( pe.EQ.from_pe )then
              if( LOC(get_data).NE.LOC(put_data) )then
!dir$ IVDEP
                  do i = 1,get_len
                     get_data(i) = put_data(i)
                  end do
              end if
          end if
          call mpp_broadcast( get_data, get_len, from_pe )
          return

      else if( to_pe.EQ.ANY_PE )then !we don't have a destination to do puts to, so only do gets

      else if( to_pe.NE.NULL_PE )then  !no other valid cases except NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid to_pe.' )
      end if

!do the get
      if( from_pe.GE.0 .AND. from_pe.LT.npes )then
          if( .NOT.allocated(local_data) ) then
             get_data_addr = LOC(get_data)
             get_len_nocomm  = get_len
          else
             do i = 1,get_len
                get_data(i) = local_data(i)
             end do
             deallocate(local_data)
          endif
      else if( from_pe.EQ.ANY_PE )then

      else if( from_pe.EQ.ALL_PES )then
          call mpp_error( FATAL, 'MPP_TRANSMIT: from_pe=ALL_PES has ambiguous meaning,' &
          & // 'and hence is not implemented.' )

      else if( from_pe.NE.NULL_PE )then !only remaining valid choice is NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid from_pe.' )
      end if

      if( debug )then
          call system_clock_default(tick)
          write( outunit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT end: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, &
                       &  put_len, get_len
      end if
      return
    end subroutine mpp_transmit_int4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_BROADCAST                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_broadcast_int4( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      integer(i4_kind), intent(inout) :: data(*)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      integer :: n

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_BROADCAST: You must first call mpp_init.' )
      return
    end subroutine mpp_broadcast_int4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_SCATTER                                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_scatterv_int4( send_data, send_counts, displs, recv_data, recv_count, root_pe, pelist, ierr)
      integer(i4_kind), dimension(:),      intent(in) :: send_data
      integer(i4_kind), dimension(:,:,:),  intent(inout) :: recv_data
      integer,                      intent(in) :: recv_count, root_pe
      integer, dimension(:),        intent(in) :: send_counts, displs, pelist
      integer,                      intent(inout) :: ierr

      integer :: n, i, j, k

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SCATTERV: You must first call mpp_init.' )

      n = 1
      do k = 1, size(recv_data, 3)
        do j = 1, size(recv_data, 2)
          do i = 1, size(recv_data, 1)
            recv_data(i,j,k) = send_data(n)
            n = n + 1
          enddo
        enddo
      enddo

    end subroutine mpp_scatterv_int4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_GATHER                                   !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_gather_int4( send_data, recv_data, count, root_pe, pelist, ierr)
      integer(i4_kind), dimension(:), intent(in) :: send_data
      integer(i4_kind), dimension(:), intent(inout) :: recv_data
      integer, dimension(:), intent(in) :: pelist(:)
      integer, intent(in) :: count, root_pe
      integer, intent(inout) :: ierr

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GATHER: You must first call mpp_init.' )

      recv_data = send_data

    end subroutine mpp_gather_int4

    subroutine mpp_gatherv_int4( send_data, send_count, recv_data, recv_counts, displs, root_pe, pelist, ierr)
      integer(i4_kind), dimension(:), intent(in) :: send_data
      integer(i4_kind), dimension(:), intent(inout) :: recv_data
      integer, dimension(:), intent(in) :: recv_counts, displs, pelist
      integer, intent(in) :: send_count, root_pe
      integer, intent(inout) :: ierr

      if( .NOT.module_is_initialized ) call mpp_error( FATAL, 'MPP_GATHERV: You must first call mpp_init.' )

      recv_data = send_data

    end subroutine mpp_gatherv_int4

!####################################################################################

# 1 "mpp/include/mpp_transmit.inc" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                  MPP_TRANSMIT                               !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_transmit_int4_scalar( put_data, to_pe, get_data, from_pe, plen, glen, block, tag, &
                                    recv_request, send_request)
      integer, intent(in) :: to_pe, from_pe
      integer(i4_kind), intent(in)  :: put_data
      integer(i4_kind), intent(out) :: get_data
      integer, optional,  intent(in) :: plen, glen
      logical, intent(in),  optional :: block
      integer, intent(in),  optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      integer                       :: put_len, get_len
      integer(i4_kind) :: put_data1D(1), get_data1D(1)
      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )

      get_data = 0

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      put_len=1; if(PRESENT(plen))put_len=plen
      get_len=1; if(PRESENT(glen))get_len=glen
      call mpp_transmit_int4 ( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                           recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_int4_scalar

    subroutine mpp_transmit_int4_2d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      integer(i4_kind), intent(in)  :: put_data(:,:)
      integer(i4_kind), intent(out) :: get_data(:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      integer(i4_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = 0

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_int4_2d

    subroutine mpp_transmit_int4_3d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      integer(i4_kind), intent(in)  :: put_data(:,:,:)
      integer(i4_kind), intent(out) :: get_data(:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      integer(i4_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = 0

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_int4_3d

    subroutine mpp_transmit_int4_4d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      integer(i4_kind), intent(in)  :: put_data(:,:,:,:)
      integer(i4_kind), intent(out) :: get_data(:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      integer(i4_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = 0

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_int4_4d

    subroutine mpp_transmit_int4_5d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      integer(i4_kind), intent(in)  :: put_data(:,:,:,:,:)
      integer(i4_kind), intent(out) :: get_data(:,:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      integer(i4_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = 0

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_int4_5d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                              MPP_SEND and RECV                              !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_recv_int4( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      integer(i4_kind), intent(out) :: get_data(*)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      integer(i4_kind) :: dummy(1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_int4

    subroutine mpp_send_int4( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      integer(i4_kind), intent(in) :: put_data(*)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      integer(i4_kind) :: dummy(1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag=tag, send_request=request )
    end subroutine mpp_send_int4

    subroutine mpp_recv_int4_scalar( get_data, from_pe, glen, block, tag, request, omp_offload )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: from_pe
      integer(i4_kind), intent(out) :: get_data
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      integer, optional, intent(in) :: glen
      logical, optional, intent(in) :: omp_offload
      integer                       :: get_len
      integer(i4_kind) :: get_data1D(1)
      integer(i4_kind) :: dummy(1)

      pointer( ptr, get_data1D )
      get_data = 0

      ptr = LOC(get_data)
      get_len=1; if(PRESENT(glen))get_len=glen
      call mpp_transmit( dummy, 1, NULL_PE, get_data1D, get_len, from_pe, &
          block, tag, recv_request=request, omp_offload=omp_offload )

    end subroutine mpp_recv_int4_scalar

    subroutine mpp_send_int4_scalar( put_data, to_pe, plen, tag, request, omp_offload)
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: to_pe
      integer(i4_kind), intent(in) :: put_data
      integer, optional, intent(in) :: plen
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      logical, optional, intent(in) :: omp_offload
      integer                       :: put_len
      integer(i4_kind) :: put_data1D(1)
      integer(i4_kind) :: dummy(1)

      pointer( ptr, put_data1D )
      ptr = LOC(put_data)
      put_len=1; if(PRESENT(plen))put_len=plen
      call mpp_transmit( put_data1D, put_len, to_pe, dummy, 1, NULL_PE, &
          tag=tag, send_request=request, omp_offload=omp_offload )

    end subroutine mpp_send_int4_scalar

    subroutine mpp_recv_int4_2d( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      integer(i4_kind), intent(out) :: get_data(:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      integer(i4_kind) :: dummy(1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, &
          block, tag, recv_request=request )
    end subroutine mpp_recv_int4_2d

    subroutine mpp_send_int4_2d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      integer(i4_kind), intent(in) :: put_data(:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      integer(i4_kind) :: dummy(1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_int4_2d

    subroutine mpp_recv_int4_3d( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      integer(i4_kind), intent(out) :: get_data(:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      integer(i4_kind) :: dummy(1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_int4_3d

    subroutine mpp_send_int4_3d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      integer(i4_kind), intent(in) :: put_data(:,:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      integer(i4_kind) :: dummy(1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_int4_3d

    subroutine mpp_recv_int4_4d( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      integer(i4_kind), intent(out) :: get_data(:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      integer(i4_kind) :: dummy(1,1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_int4_4d

    subroutine mpp_send_int4_4d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      integer(i4_kind), intent(in) :: put_data(:,:,:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      integer(i4_kind) :: dummy(1,1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_int4_4d

    subroutine mpp_recv_int4_5d( get_data, get_len, from_pe, block, tag, request)
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      integer(i4_kind), intent(out) :: get_data(:,:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      integer(i4_kind) :: dummy(1,1,1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_int4_5d

    subroutine mpp_send_int4_5d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      integer(i4_kind), intent(in) :: put_data(:,:,:,:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      integer(i4_kind) :: dummy(1,1,1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_int4_5d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_BROADCAST                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_broadcast_int4_scalar( broadcast_data, from_pe, pelist )
      integer(i4_kind), intent(inout) :: broadcast_data
      integer, intent(in) :: from_pe
      integer, intent(in), optional :: pelist(:)
      integer(i4_kind) :: data1D(1)

      pointer( ptr, data1D )

      ptr = LOC(broadcast_data)
      call mpp_broadcast_int4( data1D, 1, from_pe, pelist )

      return
    end subroutine mpp_broadcast_int4_scalar

    subroutine mpp_broadcast_int4_2d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      integer(i4_kind), intent(inout) :: broadcast_data(:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      integer(i4_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine mpp_broadcast_int4_2d

    subroutine mpp_broadcast_int4_3d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      integer(i4_kind), intent(inout) :: broadcast_data(:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      integer(i4_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
   end subroutine mpp_broadcast_int4_3d

    subroutine mpp_broadcast_int4_4d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      integer(i4_kind), intent(inout) :: broadcast_data(:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      integer(i4_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine mpp_broadcast_int4_4d

    subroutine mpp_broadcast_int4_5d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      integer(i4_kind), intent(inout) :: broadcast_data(:,:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      integer(i4_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine mpp_broadcast_int4_5d
# 217 "mpp/include/mpp_transmit_nocomm.fh" 2
# 640 "mpp/include/mpp_comm_nocomm.inc" 2

# 704 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_transmit_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                  MPP_TRANSMIT                               !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> A message-passing routine intended to be reminiscent equally of both MPI and SHMEM
!! put_data and get_data are contiguous logical(l8_kind) arrays
!!at each call, your put_data array is put to   to_pe's get_data
!!              your get_data array is got from from_pe's put_data
!!i.e we assume that typically (e.g updating halo regions) each PE performs a put _and_ a get
!!special PE designations:
!!      NULL_PE: to disable a put or a get (e.g at boundaries)
!!      ANY_PE:  if remote PE for the put or get is to be unspecific
!!      ALL_PES: broadcast and collect operations (collect not yet implemented)
!!ideally we would not pass length, but this f77-style call performs better
!!(arrays passed by address, not descriptor) further, this permits <length> contiguous
!!words from an array of any rank to be passed (avoiding f90 rank conformance check)
!!caller is responsible for completion checks (mpp_sync_self) before and after
    subroutine mpp_transmit_logical8( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, recv_request, &
                            &  send_request, omp_offload )

      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      logical(l8_kind), intent(in)  :: put_data(*)
      logical(l8_kind), intent(out) :: get_data(*)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      logical, intent(in), optional :: omp_offload
        ! NOTE: omp_offload is unused in this function

      integer :: i, outunit
      logical(l8_kind), allocatable, save :: local_data(:) !local copy used by non-parallel code (no SHMEM or MPI)
      integer(i8_kind),     save :: get_data_addr=-9999
      logical(l8_kind)                    :: get_data_local(get_len_nocomm)
      pointer(ptr_get_data, get_data_local)


      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_TRANSMIT: You must first call mpp_init.' )
      if( to_pe.EQ.NULL_PE .AND. from_pe.EQ.NULL_PE )return


      outunit = stdout()
      if( debug )then
          call system_clock_default(tick)
          write( outunit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT begin: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, &
                       &  put_len, get_len
      end if

!do put first and then get
      if( to_pe.GE.0 .AND. to_pe.LT.npes )then
          if( allocated(local_data) ) &
               call mpp_error( FATAL, 'MPP_TRANSMIT: local_data should have been deallocated by prior receive.' )
          if( get_len_nocomm > 0) then  ! pre-post recv
             ptr_get_data = get_data_addr
             do i = 1,get_len_nocomm
                get_data_local(i) = put_data(i)
             end do
             get_len_nocomm = 0
             get_data_addr = -9999
          else
             allocate( local_data(put_len) )
             do i = 1,put_len
                local_data(i) = put_data(i)
             end do
          endif
      else if( to_pe.EQ.ALL_PES )then !this is a broadcast from from_pe
          if( from_pe.LT.0 .OR. from_pe.GE.npes ) &
            call mpp_error( FATAL, 'MPP_TRANSMIT: broadcasting from invalid PE.' )
          if( put_len.GT.get_len ) &
            call mpp_error( FATAL, 'MPP_TRANSMIT: size mismatch between put_data and get_data.' )
          if( pe.EQ.from_pe )then
              if( LOC(get_data).NE.LOC(put_data) )then
!dir$ IVDEP
                  do i = 1,get_len
                     get_data(i) = put_data(i)
                  end do
              end if
          end if
          call mpp_broadcast( get_data, get_len, from_pe )
          return

      else if( to_pe.EQ.ANY_PE )then !we don't have a destination to do puts to, so only do gets

      else if( to_pe.NE.NULL_PE )then  !no other valid cases except NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid to_pe.' )
      end if

!do the get
      if( from_pe.GE.0 .AND. from_pe.LT.npes )then
          if( .NOT.allocated(local_data) ) then
             get_data_addr = LOC(get_data)
             get_len_nocomm  = get_len
          else
             do i = 1,get_len
                get_data(i) = local_data(i)
             end do
             deallocate(local_data)
          endif
      else if( from_pe.EQ.ANY_PE )then

      else if( from_pe.EQ.ALL_PES )then
          call mpp_error( FATAL, 'MPP_TRANSMIT: from_pe=ALL_PES has ambiguous meaning,' &
          & // 'and hence is not implemented.' )

      else if( from_pe.NE.NULL_PE )then !only remaining valid choice is NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid from_pe.' )
      end if

      if( debug )then
          call system_clock_default(tick)
          write( outunit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT end: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, &
                       &  put_len, get_len
      end if
      return
    end subroutine mpp_transmit_logical8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_BROADCAST                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_broadcast_logical8( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      logical(l8_kind), intent(inout) :: data(*)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      integer :: n

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_BROADCAST: You must first call mpp_init.' )
      return
    end subroutine mpp_broadcast_logical8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_SCATTER                                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_scatterv_logical8( send_data, send_counts, displs, recv_data, recv_count, root_pe, pelist, ierr)
      logical(l8_kind), dimension(:),      intent(in) :: send_data
      logical(l8_kind), dimension(:,:,:),  intent(inout) :: recv_data
      integer,                      intent(in) :: recv_count, root_pe
      integer, dimension(:),        intent(in) :: send_counts, displs, pelist
      integer,                      intent(inout) :: ierr

      integer :: n, i, j, k

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SCATTERV: You must first call mpp_init.' )

      n = 1
      do k = 1, size(recv_data, 3)
        do j = 1, size(recv_data, 2)
          do i = 1, size(recv_data, 1)
            recv_data(i,j,k) = send_data(n)
            n = n + 1
          enddo
        enddo
      enddo

    end subroutine mpp_scatterv_logical8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_GATHER                                   !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_gather_logical8( send_data, recv_data, count, root_pe, pelist, ierr)
      logical(l8_kind), dimension(:), intent(in) :: send_data
      logical(l8_kind), dimension(:), intent(inout) :: recv_data
      integer, dimension(:), intent(in) :: pelist(:)
      integer, intent(in) :: count, root_pe
      integer, intent(inout) :: ierr

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GATHER: You must first call mpp_init.' )

      recv_data = send_data

    end subroutine mpp_gather_logical8

    subroutine mpp_gatherv_logical8( send_data, send_count, recv_data, recv_counts, displs, root_pe, pelist, ierr)
      logical(l8_kind), dimension(:), intent(in) :: send_data
      logical(l8_kind), dimension(:), intent(inout) :: recv_data
      integer, dimension(:), intent(in) :: recv_counts, displs, pelist
      integer, intent(in) :: send_count, root_pe
      integer, intent(inout) :: ierr

      if( .NOT.module_is_initialized ) call mpp_error( FATAL, 'MPP_GATHERV: You must first call mpp_init.' )

      recv_data = send_data

    end subroutine mpp_gatherv_logical8

!####################################################################################

# 1 "mpp/include/mpp_transmit.inc" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                  MPP_TRANSMIT                               !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_transmit_logical8_scalar( put_data, to_pe, get_data, from_pe, plen, glen, block, tag, &
                                    recv_request, send_request)
      integer, intent(in) :: to_pe, from_pe
      logical(l8_kind), intent(in)  :: put_data
      logical(l8_kind), intent(out) :: get_data
      integer, optional,  intent(in) :: plen, glen
      logical, intent(in),  optional :: block
      integer, intent(in),  optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      integer                       :: put_len, get_len
      logical(l8_kind) :: put_data1D(1), get_data1D(1)
      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )

      get_data = .false.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      put_len=1; if(PRESENT(plen))put_len=plen
      get_len=1; if(PRESENT(glen))get_len=glen
      call mpp_transmit_logical8 ( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                           recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_logical8_scalar

    subroutine mpp_transmit_logical8_2d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      logical(l8_kind), intent(in)  :: put_data(:,:)
      logical(l8_kind), intent(out) :: get_data(:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      logical(l8_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = .false.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_logical8_2d

    subroutine mpp_transmit_logical8_3d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      logical(l8_kind), intent(in)  :: put_data(:,:,:)
      logical(l8_kind), intent(out) :: get_data(:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      logical(l8_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = .false.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_logical8_3d

    subroutine mpp_transmit_logical8_4d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      logical(l8_kind), intent(in)  :: put_data(:,:,:,:)
      logical(l8_kind), intent(out) :: get_data(:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      logical(l8_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = .false.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_logical8_4d

    subroutine mpp_transmit_logical8_5d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      logical(l8_kind), intent(in)  :: put_data(:,:,:,:,:)
      logical(l8_kind), intent(out) :: get_data(:,:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      logical(l8_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = .false.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_logical8_5d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                              MPP_SEND and RECV                              !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_recv_logical8( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      logical(l8_kind), intent(out) :: get_data(*)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      logical(l8_kind) :: dummy(1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_logical8

    subroutine mpp_send_logical8( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      logical(l8_kind), intent(in) :: put_data(*)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      logical(l8_kind) :: dummy(1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag=tag, send_request=request )
    end subroutine mpp_send_logical8

    subroutine mpp_recv_logical8_scalar( get_data, from_pe, glen, block, tag, request, omp_offload )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: from_pe
      logical(l8_kind), intent(out) :: get_data
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      integer, optional, intent(in) :: glen
      logical, optional, intent(in) :: omp_offload
      integer                       :: get_len
      logical(l8_kind) :: get_data1D(1)
      logical(l8_kind) :: dummy(1)

      pointer( ptr, get_data1D )
      get_data = .false.

      ptr = LOC(get_data)
      get_len=1; if(PRESENT(glen))get_len=glen
      call mpp_transmit( dummy, 1, NULL_PE, get_data1D, get_len, from_pe, &
          block, tag, recv_request=request, omp_offload=omp_offload )

    end subroutine mpp_recv_logical8_scalar

    subroutine mpp_send_logical8_scalar( put_data, to_pe, plen, tag, request, omp_offload)
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: to_pe
      logical(l8_kind), intent(in) :: put_data
      integer, optional, intent(in) :: plen
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      logical, optional, intent(in) :: omp_offload
      integer                       :: put_len
      logical(l8_kind) :: put_data1D(1)
      logical(l8_kind) :: dummy(1)

      pointer( ptr, put_data1D )
      ptr = LOC(put_data)
      put_len=1; if(PRESENT(plen))put_len=plen
      call mpp_transmit( put_data1D, put_len, to_pe, dummy, 1, NULL_PE, &
          tag=tag, send_request=request, omp_offload=omp_offload )

    end subroutine mpp_send_logical8_scalar

    subroutine mpp_recv_logical8_2d( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      logical(l8_kind), intent(out) :: get_data(:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      logical(l8_kind) :: dummy(1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, &
          block, tag, recv_request=request )
    end subroutine mpp_recv_logical8_2d

    subroutine mpp_send_logical8_2d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      logical(l8_kind), intent(in) :: put_data(:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      logical(l8_kind) :: dummy(1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_logical8_2d

    subroutine mpp_recv_logical8_3d( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      logical(l8_kind), intent(out) :: get_data(:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      logical(l8_kind) :: dummy(1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_logical8_3d

    subroutine mpp_send_logical8_3d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      logical(l8_kind), intent(in) :: put_data(:,:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      logical(l8_kind) :: dummy(1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_logical8_3d

    subroutine mpp_recv_logical8_4d( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      logical(l8_kind), intent(out) :: get_data(:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      logical(l8_kind) :: dummy(1,1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_logical8_4d

    subroutine mpp_send_logical8_4d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      logical(l8_kind), intent(in) :: put_data(:,:,:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      logical(l8_kind) :: dummy(1,1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_logical8_4d

    subroutine mpp_recv_logical8_5d( get_data, get_len, from_pe, block, tag, request)
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      logical(l8_kind), intent(out) :: get_data(:,:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      logical(l8_kind) :: dummy(1,1,1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_logical8_5d

    subroutine mpp_send_logical8_5d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      logical(l8_kind), intent(in) :: put_data(:,:,:,:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      logical(l8_kind) :: dummy(1,1,1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_logical8_5d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_BROADCAST                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_broadcast_logical8_scalar( broadcast_data, from_pe, pelist )
      logical(l8_kind), intent(inout) :: broadcast_data
      integer, intent(in) :: from_pe
      integer, intent(in), optional :: pelist(:)
      logical(l8_kind) :: data1D(1)

      pointer( ptr, data1D )

      ptr = LOC(broadcast_data)
      call mpp_broadcast_logical8( data1D, 1, from_pe, pelist )

      return
    end subroutine mpp_broadcast_logical8_scalar

    subroutine mpp_broadcast_logical8_2d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      logical(l8_kind), intent(inout) :: broadcast_data(:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      logical(l8_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine mpp_broadcast_logical8_2d

    subroutine mpp_broadcast_logical8_3d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      logical(l8_kind), intent(inout) :: broadcast_data(:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      logical(l8_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
   end subroutine mpp_broadcast_logical8_3d

    subroutine mpp_broadcast_logical8_4d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      logical(l8_kind), intent(inout) :: broadcast_data(:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      logical(l8_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine mpp_broadcast_logical8_4d

    subroutine mpp_broadcast_logical8_5d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      logical(l8_kind), intent(inout) :: broadcast_data(:,:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      logical(l8_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine mpp_broadcast_logical8_5d
# 217 "mpp/include/mpp_transmit_nocomm.fh" 2
# 704 "mpp/include/mpp_comm_nocomm.inc" 2

# 766 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_transmit_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                  MPP_TRANSMIT                               !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> A message-passing routine intended to be reminiscent equally of both MPI and SHMEM
!! put_data and get_data are contiguous logical(l4_kind) arrays
!!at each call, your put_data array is put to   to_pe's get_data
!!              your get_data array is got from from_pe's put_data
!!i.e we assume that typically (e.g updating halo regions) each PE performs a put _and_ a get
!!special PE designations:
!!      NULL_PE: to disable a put or a get (e.g at boundaries)
!!      ANY_PE:  if remote PE for the put or get is to be unspecific
!!      ALL_PES: broadcast and collect operations (collect not yet implemented)
!!ideally we would not pass length, but this f77-style call performs better
!!(arrays passed by address, not descriptor) further, this permits <length> contiguous
!!words from an array of any rank to be passed (avoiding f90 rank conformance check)
!!caller is responsible for completion checks (mpp_sync_self) before and after
    subroutine mpp_transmit_logical4( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, recv_request, &
                            &  send_request, omp_offload )

      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      logical(l4_kind), intent(in)  :: put_data(*)
      logical(l4_kind), intent(out) :: get_data(*)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      logical, intent(in), optional :: omp_offload
        ! NOTE: omp_offload is unused in this function

      integer :: i, outunit
      logical(l4_kind), allocatable, save :: local_data(:) !local copy used by non-parallel code (no SHMEM or MPI)
      integer(i8_kind),     save :: get_data_addr=-9999
      logical(l4_kind)                    :: get_data_local(get_len_nocomm)
      pointer(ptr_get_data, get_data_local)


      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_TRANSMIT: You must first call mpp_init.' )
      if( to_pe.EQ.NULL_PE .AND. from_pe.EQ.NULL_PE )return


      outunit = stdout()
      if( debug )then
          call system_clock_default(tick)
          write( outunit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT begin: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, &
                       &  put_len, get_len
      end if

!do put first and then get
      if( to_pe.GE.0 .AND. to_pe.LT.npes )then
          if( allocated(local_data) ) &
               call mpp_error( FATAL, 'MPP_TRANSMIT: local_data should have been deallocated by prior receive.' )
          if( get_len_nocomm > 0) then  ! pre-post recv
             ptr_get_data = get_data_addr
             do i = 1,get_len_nocomm
                get_data_local(i) = put_data(i)
             end do
             get_len_nocomm = 0
             get_data_addr = -9999
          else
             allocate( local_data(put_len) )
             do i = 1,put_len
                local_data(i) = put_data(i)
             end do
          endif
      else if( to_pe.EQ.ALL_PES )then !this is a broadcast from from_pe
          if( from_pe.LT.0 .OR. from_pe.GE.npes ) &
            call mpp_error( FATAL, 'MPP_TRANSMIT: broadcasting from invalid PE.' )
          if( put_len.GT.get_len ) &
            call mpp_error( FATAL, 'MPP_TRANSMIT: size mismatch between put_data and get_data.' )
          if( pe.EQ.from_pe )then
              if( LOC(get_data).NE.LOC(put_data) )then
!dir$ IVDEP
                  do i = 1,get_len
                     get_data(i) = put_data(i)
                  end do
              end if
          end if
          call mpp_broadcast( get_data, get_len, from_pe )
          return

      else if( to_pe.EQ.ANY_PE )then !we don't have a destination to do puts to, so only do gets

      else if( to_pe.NE.NULL_PE )then  !no other valid cases except NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid to_pe.' )
      end if

!do the get
      if( from_pe.GE.0 .AND. from_pe.LT.npes )then
          if( .NOT.allocated(local_data) ) then
             get_data_addr = LOC(get_data)
             get_len_nocomm  = get_len
          else
             do i = 1,get_len
                get_data(i) = local_data(i)
             end do
             deallocate(local_data)
          endif
      else if( from_pe.EQ.ANY_PE )then

      else if( from_pe.EQ.ALL_PES )then
          call mpp_error( FATAL, 'MPP_TRANSMIT: from_pe=ALL_PES has ambiguous meaning,' &
          & // 'and hence is not implemented.' )

      else if( from_pe.NE.NULL_PE )then !only remaining valid choice is NULL_PE
          call mpp_error( FATAL, 'MPP_TRANSMIT: invalid from_pe.' )
      end if

      if( debug )then
          call system_clock_default(tick)
          write( outunit,'(a,i18,a,i5,a,2i5,2i8)' )&
               'T=',tick, ' PE=',pe, ' MPP_TRANSMIT end: to_pe, from_pe, put_len, get_len=', to_pe, from_pe, &
                       &  put_len, get_len
      end if
      return
    end subroutine mpp_transmit_logical4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_BROADCAST                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_broadcast_logical4( data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      logical(l4_kind), intent(inout) :: data(*)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      integer :: n

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_BROADCAST: You must first call mpp_init.' )
      return
    end subroutine mpp_broadcast_logical4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_SCATTER                                  !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_scatterv_logical4( send_data, send_counts, displs, recv_data, recv_count, root_pe, pelist, ierr)
      logical(l4_kind), dimension(:),      intent(in) :: send_data
      logical(l4_kind), dimension(:,:,:),  intent(inout) :: recv_data
      integer,                      intent(in) :: recv_count, root_pe
      integer, dimension(:),        intent(in) :: send_counts, displs, pelist
      integer,                      intent(inout) :: ierr

      integer :: n, i, j, k

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SCATTERV: You must first call mpp_init.' )

      n = 1
      do k = 1, size(recv_data, 3)
        do j = 1, size(recv_data, 2)
          do i = 1, size(recv_data, 1)
            recv_data(i,j,k) = send_data(n)
            n = n + 1
          enddo
        enddo
      enddo

    end subroutine mpp_scatterv_logical4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_GATHER                                   !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mpp_gather_logical4( send_data, recv_data, count, root_pe, pelist, ierr)
      logical(l4_kind), dimension(:), intent(in) :: send_data
      logical(l4_kind), dimension(:), intent(inout) :: recv_data
      integer, dimension(:), intent(in) :: pelist(:)
      integer, intent(in) :: count, root_pe
      integer, intent(inout) :: ierr

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GATHER: You must first call mpp_init.' )

      recv_data = send_data

    end subroutine mpp_gather_logical4

    subroutine mpp_gatherv_logical4( send_data, send_count, recv_data, recv_counts, displs, root_pe, pelist, ierr)
      logical(l4_kind), dimension(:), intent(in) :: send_data
      logical(l4_kind), dimension(:), intent(inout) :: recv_data
      integer, dimension(:), intent(in) :: recv_counts, displs, pelist
      integer, intent(in) :: send_count, root_pe
      integer, intent(inout) :: ierr

      if( .NOT.module_is_initialized ) call mpp_error( FATAL, 'MPP_GATHERV: You must first call mpp_init.' )

      recv_data = send_data

    end subroutine mpp_gatherv_logical4

!####################################################################################

# 1 "mpp/include/mpp_transmit.inc" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                  MPP_TRANSMIT                               !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_transmit_logical4_scalar( put_data, to_pe, get_data, from_pe, plen, glen, block, tag, &
                                    recv_request, send_request)
      integer, intent(in) :: to_pe, from_pe
      logical(l4_kind), intent(in)  :: put_data
      logical(l4_kind), intent(out) :: get_data
      integer, optional,  intent(in) :: plen, glen
      logical, intent(in),  optional :: block
      integer, intent(in),  optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      integer                       :: put_len, get_len
      logical(l4_kind) :: put_data1D(1), get_data1D(1)
      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )

      get_data = .false.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      put_len=1; if(PRESENT(plen))put_len=plen
      get_len=1; if(PRESENT(glen))get_len=glen
      call mpp_transmit_logical4 ( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                           recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_logical4_scalar

    subroutine mpp_transmit_logical4_2d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      logical(l4_kind), intent(in)  :: put_data(:,:)
      logical(l4_kind), intent(out) :: get_data(:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      logical(l4_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = .false.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_logical4_2d

    subroutine mpp_transmit_logical4_3d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      logical(l4_kind), intent(in)  :: put_data(:,:,:)
      logical(l4_kind), intent(out) :: get_data(:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      logical(l4_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = .false.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_logical4_3d

    subroutine mpp_transmit_logical4_4d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      logical(l4_kind), intent(in)  :: put_data(:,:,:,:)
      logical(l4_kind), intent(out) :: get_data(:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      logical(l4_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = .false.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_logical4_4d

    subroutine mpp_transmit_logical4_5d( put_data, put_len, to_pe, get_data, get_len, from_pe, block, tag, &
                                recv_request, send_request )
      integer, intent(in) :: put_len, to_pe, get_len, from_pe
      logical(l4_kind), intent(in)  :: put_data(:,:,:,:,:)
      logical(l4_kind), intent(out) :: get_data(:,:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: recv_request, send_request
      logical(l4_kind) :: put_data1D(put_len), get_data1D(get_len)

      pointer( ptrp, put_data1D )
      pointer( ptrg, get_data1D )
      get_data = .false.

      ptrp = LOC(put_data)
      ptrg = LOC(get_data)
      call mpp_transmit( put_data1D, put_len, to_pe, get_data1D, get_len, from_pe, block, tag, &
                         recv_request=recv_request, send_request=send_request )

      return
    end subroutine mpp_transmit_logical4_5d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                              MPP_SEND and RECV                              !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_recv_logical4( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      logical(l4_kind), intent(out) :: get_data(*)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      logical(l4_kind) :: dummy(1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_logical4

    subroutine mpp_send_logical4( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      logical(l4_kind), intent(in) :: put_data(*)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      logical(l4_kind) :: dummy(1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag=tag, send_request=request )
    end subroutine mpp_send_logical4

    subroutine mpp_recv_logical4_scalar( get_data, from_pe, glen, block, tag, request, omp_offload )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: from_pe
      logical(l4_kind), intent(out) :: get_data
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      integer, optional, intent(in) :: glen
      logical, optional, intent(in) :: omp_offload
      integer                       :: get_len
      logical(l4_kind) :: get_data1D(1)
      logical(l4_kind) :: dummy(1)

      pointer( ptr, get_data1D )
      get_data = .false.

      ptr = LOC(get_data)
      get_len=1; if(PRESENT(glen))get_len=glen
      call mpp_transmit( dummy, 1, NULL_PE, get_data1D, get_len, from_pe, &
          block, tag, recv_request=request, omp_offload=omp_offload )

    end subroutine mpp_recv_logical4_scalar

    subroutine mpp_send_logical4_scalar( put_data, to_pe, plen, tag, request, omp_offload)
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: to_pe
      logical(l4_kind), intent(in) :: put_data
      integer, optional, intent(in) :: plen
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      logical, optional, intent(in) :: omp_offload
      integer                       :: put_len
      logical(l4_kind) :: put_data1D(1)
      logical(l4_kind) :: dummy(1)

      pointer( ptr, put_data1D )
      ptr = LOC(put_data)
      put_len=1; if(PRESENT(plen))put_len=plen
      call mpp_transmit( put_data1D, put_len, to_pe, dummy, 1, NULL_PE, &
          tag=tag, send_request=request, omp_offload=omp_offload )

    end subroutine mpp_send_logical4_scalar

    subroutine mpp_recv_logical4_2d( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      logical(l4_kind), intent(out) :: get_data(:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      logical(l4_kind) :: dummy(1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, &
          block, tag, recv_request=request )
    end subroutine mpp_recv_logical4_2d

    subroutine mpp_send_logical4_2d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      logical(l4_kind), intent(in) :: put_data(:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      logical(l4_kind) :: dummy(1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_logical4_2d

    subroutine mpp_recv_logical4_3d( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      logical(l4_kind), intent(out) :: get_data(:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      logical(l4_kind) :: dummy(1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_logical4_3d

    subroutine mpp_send_logical4_3d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      logical(l4_kind), intent(in) :: put_data(:,:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      logical(l4_kind) :: dummy(1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_logical4_3d

    subroutine mpp_recv_logical4_4d( get_data, get_len, from_pe, block, tag, request )
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      logical(l4_kind), intent(out) :: get_data(:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      logical(l4_kind) :: dummy(1,1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_logical4_4d

    subroutine mpp_send_logical4_4d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      logical(l4_kind), intent(in) :: put_data(:,:,:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      logical(l4_kind) :: dummy(1,1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_logical4_4d

    subroutine mpp_recv_logical4_5d( get_data, get_len, from_pe, block, tag, request)
!a mpp_transmit with null arguments on the put side
      integer, intent(in) :: get_len, from_pe
      logical(l4_kind), intent(out) :: get_data(:,:,:,:,:)
      logical, intent(in), optional :: block
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request

      logical(l4_kind) :: dummy(1,1,1,1,1)
      call mpp_transmit( dummy, 1, NULL_PE, get_data, get_len, from_pe, block, tag, recv_request=request )
    end subroutine mpp_recv_logical4_5d

    subroutine mpp_send_logical4_5d( put_data, put_len, to_pe, tag, request )
!a mpp_transmit with null arguments on the get side
      integer, intent(in) :: put_len, to_pe
      logical(l4_kind), intent(in) :: put_data(:,:,:,:,:)
      integer, intent(in), optional :: tag
      type(mpi_request), intent(out), optional :: request
      logical(l4_kind) :: dummy(1,1,1,1,1)
      call mpp_transmit( put_data, put_len, to_pe, dummy, 1, NULL_PE, tag = tag, send_request=request )
    end subroutine mpp_send_logical4_5d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_BROADCAST                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpp_broadcast_logical4_scalar( broadcast_data, from_pe, pelist )
      logical(l4_kind), intent(inout) :: broadcast_data
      integer, intent(in) :: from_pe
      integer, intent(in), optional :: pelist(:)
      logical(l4_kind) :: data1D(1)

      pointer( ptr, data1D )

      ptr = LOC(broadcast_data)
      call mpp_broadcast_logical4( data1D, 1, from_pe, pelist )

      return
    end subroutine mpp_broadcast_logical4_scalar

    subroutine mpp_broadcast_logical4_2d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      logical(l4_kind), intent(inout) :: broadcast_data(:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      logical(l4_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine mpp_broadcast_logical4_2d

    subroutine mpp_broadcast_logical4_3d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      logical(l4_kind), intent(inout) :: broadcast_data(:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      logical(l4_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
   end subroutine mpp_broadcast_logical4_3d

    subroutine mpp_broadcast_logical4_4d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      logical(l4_kind), intent(inout) :: broadcast_data(:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      logical(l4_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine mpp_broadcast_logical4_4d

    subroutine mpp_broadcast_logical4_5d( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      logical(l4_kind), intent(inout) :: broadcast_data(:,:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      logical(l4_kind) :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine mpp_broadcast_logical4_5d
# 217 "mpp/include/mpp_transmit_nocomm.fh" 2
# 766 "mpp/include/mpp_comm_nocomm.inc" 2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!            GLOBAL REDUCTION ROUTINES: mpp_max, mpp_sum, mpp_min             !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 786 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_reduce_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
    !> Find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast to all PEs. Nocomm version
    subroutine mpp_max_real8_0d( a, pelist )
      real(r8_kind), intent(inout) :: a
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_0D: You must first call mpp_init.' )
      return
    end subroutine mpp_max_real8_0d

    !> Find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast to all PEs. Nocomm version
    subroutine mpp_max_real8_1d( a, length, pelist )
      real(r8_kind), intent(inout) :: a(:)
      integer,   intent(in)    :: length
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_1D: You must first call mpp_init.' )
      return
    end subroutine mpp_max_real8_1d

# 786 "mpp/include/mpp_comm_nocomm.inc" 2

# 800 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_reduce_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
    !> Find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast to all PEs. Nocomm version
    subroutine mpp_max_real4_0d( a, pelist )
      real(r4_kind), intent(inout) :: a
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_0D: You must first call mpp_init.' )
      return
    end subroutine mpp_max_real4_0d

    !> Find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast to all PEs. Nocomm version
    subroutine mpp_max_real4_1d( a, length, pelist )
      real(r4_kind), intent(inout) :: a(:)
      integer,   intent(in)    :: length
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_1D: You must first call mpp_init.' )
      return
    end subroutine mpp_max_real4_1d

# 800 "mpp/include/mpp_comm_nocomm.inc" 2

# 814 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_reduce_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
    !> Find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast to all PEs. Nocomm version
    subroutine mpp_max_int8_0d( a, pelist )
      integer(i8_kind), intent(inout) :: a
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_0D: You must first call mpp_init.' )
      return
    end subroutine mpp_max_int8_0d

    !> Find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast to all PEs. Nocomm version
    subroutine mpp_max_int8_1d( a, length, pelist )
      integer(i8_kind), intent(inout) :: a(:)
      integer,   intent(in)    :: length
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_1D: You must first call mpp_init.' )
      return
    end subroutine mpp_max_int8_1d

# 814 "mpp/include/mpp_comm_nocomm.inc" 2

# 828 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_reduce_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
    !> Find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast to all PEs. Nocomm version
    subroutine mpp_max_int4_0d( a, pelist )
      integer(i4_kind), intent(inout) :: a
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_0D: You must first call mpp_init.' )
      return
    end subroutine mpp_max_int4_0d

    !> Find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast to all PEs. Nocomm version
    subroutine mpp_max_int4_1d( a, length, pelist )
      integer(i4_kind), intent(inout) :: a(:)
      integer,   intent(in)    :: length
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_1D: You must first call mpp_init.' )
      return
    end subroutine mpp_max_int4_1d

# 828 "mpp/include/mpp_comm_nocomm.inc" 2

# 842 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_reduce_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
    !> Find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast to all PEs. Nocomm version
    subroutine mpp_min_real8_0d( a, pelist )
      real(r8_kind), intent(inout) :: a
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_0D: You must first call mpp_init.' )
      return
    end subroutine mpp_min_real8_0d

    !> Find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast to all PEs. Nocomm version
    subroutine mpp_min_real8_1d( a, length, pelist )
      real(r8_kind), intent(inout) :: a(:)
      integer,   intent(in)    :: length
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_1D: You must first call mpp_init.' )
      return
    end subroutine mpp_min_real8_1d

# 842 "mpp/include/mpp_comm_nocomm.inc" 2

# 856 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_reduce_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
    !> Find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast to all PEs. Nocomm version
    subroutine mpp_min_real4_0d( a, pelist )
      real(r4_kind), intent(inout) :: a
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_0D: You must first call mpp_init.' )
      return
    end subroutine mpp_min_real4_0d

    !> Find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast to all PEs. Nocomm version
    subroutine mpp_min_real4_1d( a, length, pelist )
      real(r4_kind), intent(inout) :: a(:)
      integer,   intent(in)    :: length
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_1D: You must first call mpp_init.' )
      return
    end subroutine mpp_min_real4_1d

# 856 "mpp/include/mpp_comm_nocomm.inc" 2

# 870 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_reduce_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
    !> Find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast to all PEs. Nocomm version
    subroutine mpp_min_int8_0d( a, pelist )
      integer(i8_kind), intent(inout) :: a
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_0D: You must first call mpp_init.' )
      return
    end subroutine mpp_min_int8_0d

    !> Find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast to all PEs. Nocomm version
    subroutine mpp_min_int8_1d( a, length, pelist )
      integer(i8_kind), intent(inout) :: a(:)
      integer,   intent(in)    :: length
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_1D: You must first call mpp_init.' )
      return
    end subroutine mpp_min_int8_1d

# 870 "mpp/include/mpp_comm_nocomm.inc" 2

# 884 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_reduce_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
    !> Find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast to all PEs. Nocomm version
    subroutine mpp_min_int4_0d( a, pelist )
      integer(i4_kind), intent(inout) :: a
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_0D: You must first call mpp_init.' )
      return
    end subroutine mpp_min_int4_0d

    !> Find the max of scalar a the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast to all PEs. Nocomm version
    subroutine mpp_min_int4_1d( a, length, pelist )
      integer(i4_kind), intent(inout) :: a(:)
      integer,   intent(in)    :: length
      integer, intent(in), optional :: pelist(0:)
      integer :: n
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_REDUCE_1D: You must first call mpp_init.' )
      return
    end subroutine mpp_min_int4_1d

# 884 "mpp/include/mpp_comm_nocomm.inc" 2

# 904 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_sum_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
    !> Sums array a over the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast: all PEs have the sum in a at the end
    !! we are using f77-style call: array passed by address and not descriptor; further,
    !! the f90 conformance check is avoided.
    subroutine mpp_sum_real8( a, length, pelist )
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      real(r8_kind), intent(inout) :: a(*)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SUM: You must first call mpp_init.' )
      return
    end subroutine mpp_sum_real8

!#######################################################################

# 1 "mpp/include/mpp_sum.inc" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!#######################################################################

    !> Sums array a when only first element is passed: this routine just converts to a call to mpp_sum_real8
    subroutine mpp_sum_real8_scalar( a, pelist )
      real(r8_kind), intent(inout) :: a
      integer, intent(in), optional :: pelist(:)
      real(r8_kind) :: b(1)

      b(1) = a
      if( debug )call mpp_error( NOTE, 'MPP_SUM_SCALAR_: calling MPP_SUM_ ...' )
      call mpp_sum_real8( b, 1, pelist )
      a = b(1)
      return
    end subroutine mpp_sum_real8_scalar

!#######################################################################
    !> Sums 2d array across pes
    subroutine mpp_sum_real8_2d( a, length, pelist )
      real(r8_kind), intent(inout) :: a(:,:) !< 2d array to sum
      integer, intent(in) :: length !< amount of indices in given 2d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      real(r8_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine mpp_sum_real8_2d

!#######################################################################
    !> Sums 3d array across pes
    subroutine mpp_sum_real8_3d( a, length, pelist )
      real(r8_kind), intent(inout) :: a(:,:,:) !< 3d array to sum
      integer, intent(in) :: length !< amount of indices in given 3d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      real(r8_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine mpp_sum_real8_3d

!#######################################################################
    !> Sums 4d array across pes
    subroutine mpp_sum_real8_4d( a, length, pelist )
      real(r8_kind), intent(inout) :: a(:,:,:,:) !< 4d array to sum
      integer, intent(in) :: length !< amount of indices in given 4d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      real(r8_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine mpp_sum_real8_4d

!#######################################################################
    !> Sums 5d array across pes
    subroutine mpp_sum_real8_5d( a, length, pelist )
      real(r8_kind), intent(inout) :: a(:,:,:,:,:) !< 5d array to sum
      integer, intent(in) :: length !< amount of indices in given 5d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      real(r8_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine mpp_sum_real8_5d
# 33 "mpp/include/mpp_sum_nocomm.fh" 2

# 904 "mpp/include/mpp_comm_nocomm.inc" 2

# 926 "mpp/include/mpp_comm_nocomm.inc"

# 946 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_sum_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
    !> Sums array a over the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast: all PEs have the sum in a at the end
    !! we are using f77-style call: array passed by address and not descriptor; further,
    !! the f90 conformance check is avoided.
    subroutine mpp_sum_real4( a, length, pelist )
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      real(r4_kind), intent(inout) :: a(*)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SUM: You must first call mpp_init.' )
      return
    end subroutine mpp_sum_real4

!#######################################################################

# 1 "mpp/include/mpp_sum.inc" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!#######################################################################

    !> Sums array a when only first element is passed: this routine just converts to a call to mpp_sum_real4
    subroutine mpp_sum_real4_scalar( a, pelist )
      real(r4_kind), intent(inout) :: a
      integer, intent(in), optional :: pelist(:)
      real(r4_kind) :: b(1)

      b(1) = a
      if( debug )call mpp_error( NOTE, 'MPP_SUM_SCALAR_: calling MPP_SUM_ ...' )
      call mpp_sum_real4( b, 1, pelist )
      a = b(1)
      return
    end subroutine mpp_sum_real4_scalar

!#######################################################################
    !> Sums 2d array across pes
    subroutine mpp_sum_real4_2d( a, length, pelist )
      real(r4_kind), intent(inout) :: a(:,:) !< 2d array to sum
      integer, intent(in) :: length !< amount of indices in given 2d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      real(r4_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine mpp_sum_real4_2d

!#######################################################################
    !> Sums 3d array across pes
    subroutine mpp_sum_real4_3d( a, length, pelist )
      real(r4_kind), intent(inout) :: a(:,:,:) !< 3d array to sum
      integer, intent(in) :: length !< amount of indices in given 3d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      real(r4_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine mpp_sum_real4_3d

!#######################################################################
    !> Sums 4d array across pes
    subroutine mpp_sum_real4_4d( a, length, pelist )
      real(r4_kind), intent(inout) :: a(:,:,:,:) !< 4d array to sum
      integer, intent(in) :: length !< amount of indices in given 4d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      real(r4_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine mpp_sum_real4_4d

!#######################################################################
    !> Sums 5d array across pes
    subroutine mpp_sum_real4_5d( a, length, pelist )
      real(r4_kind), intent(inout) :: a(:,:,:,:,:) !< 5d array to sum
      integer, intent(in) :: length !< amount of indices in given 5d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      real(r4_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine mpp_sum_real4_5d
# 33 "mpp/include/mpp_sum_nocomm.fh" 2

# 946 "mpp/include/mpp_comm_nocomm.inc" 2

# 968 "mpp/include/mpp_comm_nocomm.inc"

# 988 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_sum_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
    !> Sums array a over the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast: all PEs have the sum in a at the end
    !! we are using f77-style call: array passed by address and not descriptor; further,
    !! the f90 conformance check is avoided.
    subroutine mpp_sum_int8( a, length, pelist )
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      integer(i8_kind), intent(inout) :: a(*)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SUM: You must first call mpp_init.' )
      return
    end subroutine mpp_sum_int8

!#######################################################################

# 1 "mpp/include/mpp_sum.inc" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!#######################################################################

    !> Sums array a when only first element is passed: this routine just converts to a call to mpp_sum_int8
    subroutine mpp_sum_int8_scalar( a, pelist )
      integer(i8_kind), intent(inout) :: a
      integer, intent(in), optional :: pelist(:)
      integer(i8_kind) :: b(1)

      b(1) = a
      if( debug )call mpp_error( NOTE, 'MPP_SUM_SCALAR_: calling MPP_SUM_ ...' )
      call mpp_sum_int8( b, 1, pelist )
      a = b(1)
      return
    end subroutine mpp_sum_int8_scalar

!#######################################################################
    !> Sums 2d array across pes
    subroutine mpp_sum_int8_2d( a, length, pelist )
      integer(i8_kind), intent(inout) :: a(:,:) !< 2d array to sum
      integer, intent(in) :: length !< amount of indices in given 2d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      integer(i8_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine mpp_sum_int8_2d

!#######################################################################
    !> Sums 3d array across pes
    subroutine mpp_sum_int8_3d( a, length, pelist )
      integer(i8_kind), intent(inout) :: a(:,:,:) !< 3d array to sum
      integer, intent(in) :: length !< amount of indices in given 3d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      integer(i8_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine mpp_sum_int8_3d

!#######################################################################
    !> Sums 4d array across pes
    subroutine mpp_sum_int8_4d( a, length, pelist )
      integer(i8_kind), intent(inout) :: a(:,:,:,:) !< 4d array to sum
      integer, intent(in) :: length !< amount of indices in given 4d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      integer(i8_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine mpp_sum_int8_4d

!#######################################################################
    !> Sums 5d array across pes
    subroutine mpp_sum_int8_5d( a, length, pelist )
      integer(i8_kind), intent(inout) :: a(:,:,:,:,:) !< 5d array to sum
      integer, intent(in) :: length !< amount of indices in given 5d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      integer(i8_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine mpp_sum_int8_5d
# 33 "mpp/include/mpp_sum_nocomm.fh" 2

# 988 "mpp/include/mpp_comm_nocomm.inc" 2

# 1008 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_sum_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
    !> Sums array a over the PEs in pelist (all PEs if this argument is omitted)
    !! result is also automatically broadcast: all PEs have the sum in a at the end
    !! we are using f77-style call: array passed by address and not descriptor; further,
    !! the f90 conformance check is avoided.
    subroutine mpp_sum_int4( a, length, pelist )
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      integer(i4_kind), intent(inout) :: a(*)

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SUM: You must first call mpp_init.' )
      return
    end subroutine mpp_sum_int4

!#######################################################################

# 1 "mpp/include/mpp_sum.inc" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!#######################################################################

    !> Sums array a when only first element is passed: this routine just converts to a call to mpp_sum_int4
    subroutine mpp_sum_int4_scalar( a, pelist )
      integer(i4_kind), intent(inout) :: a
      integer, intent(in), optional :: pelist(:)
      integer(i4_kind) :: b(1)

      b(1) = a
      if( debug )call mpp_error( NOTE, 'MPP_SUM_SCALAR_: calling MPP_SUM_ ...' )
      call mpp_sum_int4( b, 1, pelist )
      a = b(1)
      return
    end subroutine mpp_sum_int4_scalar

!#######################################################################
    !> Sums 2d array across pes
    subroutine mpp_sum_int4_2d( a, length, pelist )
      integer(i4_kind), intent(inout) :: a(:,:) !< 2d array to sum
      integer, intent(in) :: length !< amount of indices in given 2d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      integer(i4_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine mpp_sum_int4_2d

!#######################################################################
    !> Sums 3d array across pes
    subroutine mpp_sum_int4_3d( a, length, pelist )
      integer(i4_kind), intent(inout) :: a(:,:,:) !< 3d array to sum
      integer, intent(in) :: length !< amount of indices in given 3d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      integer(i4_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine mpp_sum_int4_3d

!#######################################################################
    !> Sums 4d array across pes
    subroutine mpp_sum_int4_4d( a, length, pelist )
      integer(i4_kind), intent(inout) :: a(:,:,:,:) !< 4d array to sum
      integer, intent(in) :: length !< amount of indices in given 4d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      integer(i4_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine mpp_sum_int4_4d

!#######################################################################
    !> Sums 5d array across pes
    subroutine mpp_sum_int4_5d( a, length, pelist )
      integer(i4_kind), intent(inout) :: a(:,:,:,:,:) !< 5d array to sum
      integer, intent(in) :: length !< amount of indices in given 5d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      integer(i4_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum( a1D, length, pelist )

      return
    end subroutine mpp_sum_int4_5d
# 33 "mpp/include/mpp_sum_nocomm.fh" 2

# 1008 "mpp/include/mpp_comm_nocomm.inc" 2
!--------------------------------
# 1028 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_sum_nocomm_ad.fh" 1
! -*-f90-*-


!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

    !> Sums array a over the PEs in pelist (all PEs if this argument is omitted)
    !! forward array is already summed and broadcasted: all PEs already have the ad sum
    !! we are using f77-style call: array passed by address and not descriptor; further,
    !! the f90 conformance check is avoided.
    subroutine mpp_sum_real8_ad( a, length, pelist )
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      real(r8_kind), intent(inout) :: a(*)
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SUM: You must first call mpp_init.' )
      return
    end subroutine mpp_sum_real8_ad

!#######################################################################

# 1 "mpp/include/mpp_sum_ad.inc" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!#######################################################################

    !> Sums array a.
    !! when only first element is passed: this routine just converts to a call to mpp_sum_int4
    subroutine mpp_sum_real8_scalar_ad( a, pelist )
      real(r8_kind), intent(inout) :: a !< scalar value to sum over pelist
      integer, intent(in), optional :: pelist(:)
      real(r8_kind) :: b(1)

      b(1) = a
      if( debug )call mpp_error( NOTE, 'MPP_SUM_SCALAR_: calling MPP_SUM_ ...' )
      call mpp_sum_real8_ad( b, 1, pelist )
      a = b(1)
      return
    end subroutine mpp_sum_real8_scalar_ad

!#######################################################################
    !> Sums 2d array across pes
    subroutine mpp_sum_real8_2d_ad( a, length, pelist )
      real(r8_kind), intent(inout) :: a(:,:) !< 2d array to sum
      integer, intent(in) :: length !< amount of indices in given 2d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      real(r8_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine mpp_sum_real8_2d_ad

!#######################################################################
    !> Sums 3d array across pes
    subroutine mpp_sum_real8_3d_ad( a, length, pelist )
      real(r8_kind), intent(inout) :: a(:,:,:) !< 3d array to sum
      integer, intent(in) :: length !< amount of indices in given 3d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      real(r8_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine mpp_sum_real8_3d_ad

!#######################################################################
    !> Sums 4d array across pes
    subroutine mpp_sum_real8_4d_ad( a, length, pelist )
      real(r8_kind), intent(inout) :: a(:,:,:,:) !< 4d array to sum
      integer, intent(in) :: length !< amount of indices in given 4d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      real(r8_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine mpp_sum_real8_4d_ad

!#######################################################################
    !> Sums 5d array across pes
    subroutine mpp_sum_real8_5d_ad( a, length, pelist )
      real(r8_kind), intent(inout) :: a(:,:,:,:,:) !< 5d array to sum
      integer, intent(in) :: length !< amount of indices in given 5d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      real(r8_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine mpp_sum_real8_5d_ad
# 36 "mpp/include/mpp_sum_nocomm_ad.fh" 2
# 1028 "mpp/include/mpp_comm_nocomm.inc" 2

# 1050 "mpp/include/mpp_comm_nocomm.inc"

# 1070 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_sum_nocomm_ad.fh" 1
! -*-f90-*-


!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

    !> Sums array a over the PEs in pelist (all PEs if this argument is omitted)
    !! forward array is already summed and broadcasted: all PEs already have the ad sum
    !! we are using f77-style call: array passed by address and not descriptor; further,
    !! the f90 conformance check is avoided.
    subroutine mpp_sum_real4_ad( a, length, pelist )
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      real(r4_kind), intent(inout) :: a(*)
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SUM: You must first call mpp_init.' )
      return
    end subroutine mpp_sum_real4_ad

!#######################################################################

# 1 "mpp/include/mpp_sum_ad.inc" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!#######################################################################

    !> Sums array a.
    !! when only first element is passed: this routine just converts to a call to mpp_sum_int4
    subroutine mpp_sum_real4_scalar_ad( a, pelist )
      real(r4_kind), intent(inout) :: a !< scalar value to sum over pelist
      integer, intent(in), optional :: pelist(:)
      real(r4_kind) :: b(1)

      b(1) = a
      if( debug )call mpp_error( NOTE, 'MPP_SUM_SCALAR_: calling MPP_SUM_ ...' )
      call mpp_sum_real4_ad( b, 1, pelist )
      a = b(1)
      return
    end subroutine mpp_sum_real4_scalar_ad

!#######################################################################
    !> Sums 2d array across pes
    subroutine mpp_sum_real4_2d_ad( a, length, pelist )
      real(r4_kind), intent(inout) :: a(:,:) !< 2d array to sum
      integer, intent(in) :: length !< amount of indices in given 2d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      real(r4_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine mpp_sum_real4_2d_ad

!#######################################################################
    !> Sums 3d array across pes
    subroutine mpp_sum_real4_3d_ad( a, length, pelist )
      real(r4_kind), intent(inout) :: a(:,:,:) !< 3d array to sum
      integer, intent(in) :: length !< amount of indices in given 3d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      real(r4_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine mpp_sum_real4_3d_ad

!#######################################################################
    !> Sums 4d array across pes
    subroutine mpp_sum_real4_4d_ad( a, length, pelist )
      real(r4_kind), intent(inout) :: a(:,:,:,:) !< 4d array to sum
      integer, intent(in) :: length !< amount of indices in given 4d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      real(r4_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine mpp_sum_real4_4d_ad

!#######################################################################
    !> Sums 5d array across pes
    subroutine mpp_sum_real4_5d_ad( a, length, pelist )
      real(r4_kind), intent(inout) :: a(:,:,:,:,:) !< 5d array to sum
      integer, intent(in) :: length !< amount of indices in given 5d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      real(r4_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine mpp_sum_real4_5d_ad
# 36 "mpp/include/mpp_sum_nocomm_ad.fh" 2
# 1070 "mpp/include/mpp_comm_nocomm.inc" 2

# 1092 "mpp/include/mpp_comm_nocomm.inc"

# 1112 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_sum_nocomm_ad.fh" 1
! -*-f90-*-


!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

    !> Sums array a over the PEs in pelist (all PEs if this argument is omitted)
    !! forward array is already summed and broadcasted: all PEs already have the ad sum
    !! we are using f77-style call: array passed by address and not descriptor; further,
    !! the f90 conformance check is avoided.
    subroutine mpp_sum_int8_ad( a, length, pelist )
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      integer(i8_kind), intent(inout) :: a(*)
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SUM: You must first call mpp_init.' )
      return
    end subroutine mpp_sum_int8_ad

!#######################################################################

# 1 "mpp/include/mpp_sum_ad.inc" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!#######################################################################

    !> Sums array a.
    !! when only first element is passed: this routine just converts to a call to mpp_sum_int4
    subroutine mpp_sum_int8_scalar_ad( a, pelist )
      integer(i8_kind), intent(inout) :: a !< scalar value to sum over pelist
      integer, intent(in), optional :: pelist(:)
      integer(i8_kind) :: b(1)

      b(1) = a
      if( debug )call mpp_error( NOTE, 'MPP_SUM_SCALAR_: calling MPP_SUM_ ...' )
      call mpp_sum_int8_ad( b, 1, pelist )
      a = b(1)
      return
    end subroutine mpp_sum_int8_scalar_ad

!#######################################################################
    !> Sums 2d array across pes
    subroutine mpp_sum_int8_2d_ad( a, length, pelist )
      integer(i8_kind), intent(inout) :: a(:,:) !< 2d array to sum
      integer, intent(in) :: length !< amount of indices in given 2d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      integer(i8_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine mpp_sum_int8_2d_ad

!#######################################################################
    !> Sums 3d array across pes
    subroutine mpp_sum_int8_3d_ad( a, length, pelist )
      integer(i8_kind), intent(inout) :: a(:,:,:) !< 3d array to sum
      integer, intent(in) :: length !< amount of indices in given 3d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      integer(i8_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine mpp_sum_int8_3d_ad

!#######################################################################
    !> Sums 4d array across pes
    subroutine mpp_sum_int8_4d_ad( a, length, pelist )
      integer(i8_kind), intent(inout) :: a(:,:,:,:) !< 4d array to sum
      integer, intent(in) :: length !< amount of indices in given 4d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      integer(i8_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine mpp_sum_int8_4d_ad

!#######################################################################
    !> Sums 5d array across pes
    subroutine mpp_sum_int8_5d_ad( a, length, pelist )
      integer(i8_kind), intent(inout) :: a(:,:,:,:,:) !< 5d array to sum
      integer, intent(in) :: length !< amount of indices in given 5d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      integer(i8_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine mpp_sum_int8_5d_ad
# 36 "mpp/include/mpp_sum_nocomm_ad.fh" 2
# 1112 "mpp/include/mpp_comm_nocomm.inc" 2

# 1132 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_sum_nocomm_ad.fh" 1
! -*-f90-*-


!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

    !> Sums array a over the PEs in pelist (all PEs if this argument is omitted)
    !! forward array is already summed and broadcasted: all PEs already have the ad sum
    !! we are using f77-style call: array passed by address and not descriptor; further,
    !! the f90 conformance check is avoided.
    subroutine mpp_sum_int4_ad( a, length, pelist )
      integer, intent(in) :: length
      integer, intent(in), optional :: pelist(:)
      integer(i4_kind), intent(inout) :: a(*)
      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_SUM: You must first call mpp_init.' )
      return
    end subroutine mpp_sum_int4_ad

!#######################################################################

# 1 "mpp/include/mpp_sum_ad.inc" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!#######################################################################

    !> Sums array a.
    !! when only first element is passed: this routine just converts to a call to mpp_sum_int4
    subroutine mpp_sum_int4_scalar_ad( a, pelist )
      integer(i4_kind), intent(inout) :: a !< scalar value to sum over pelist
      integer, intent(in), optional :: pelist(:)
      integer(i4_kind) :: b(1)

      b(1) = a
      if( debug )call mpp_error( NOTE, 'MPP_SUM_SCALAR_: calling MPP_SUM_ ...' )
      call mpp_sum_int4_ad( b, 1, pelist )
      a = b(1)
      return
    end subroutine mpp_sum_int4_scalar_ad

!#######################################################################
    !> Sums 2d array across pes
    subroutine mpp_sum_int4_2d_ad( a, length, pelist )
      integer(i4_kind), intent(inout) :: a(:,:) !< 2d array to sum
      integer, intent(in) :: length !< amount of indices in given 2d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      integer(i4_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine mpp_sum_int4_2d_ad

!#######################################################################
    !> Sums 3d array across pes
    subroutine mpp_sum_int4_3d_ad( a, length, pelist )
      integer(i4_kind), intent(inout) :: a(:,:,:) !< 3d array to sum
      integer, intent(in) :: length !< amount of indices in given 3d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      integer(i4_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine mpp_sum_int4_3d_ad

!#######################################################################
    !> Sums 4d array across pes
    subroutine mpp_sum_int4_4d_ad( a, length, pelist )
      integer(i4_kind), intent(inout) :: a(:,:,:,:) !< 4d array to sum
      integer, intent(in) :: length !< amount of indices in given 4d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      integer(i4_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine mpp_sum_int4_4d_ad

!#######################################################################
    !> Sums 5d array across pes
    subroutine mpp_sum_int4_5d_ad( a, length, pelist )
      integer(i4_kind), intent(inout) :: a(:,:,:,:,:) !< 5d array to sum
      integer, intent(in) :: length !< amount of indices in given 5d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      integer(i4_kind) :: a1D(length)

      pointer( ptr, a1D )
      ptr = LOC(a)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine mpp_sum_int4_5d_ad
# 36 "mpp/include/mpp_sum_nocomm_ad.fh" 2
# 1132 "mpp/include/mpp_comm_nocomm.inc" 2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!            SCATTER AND GATHER ROUTINES: mpp_alltoall                        !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# 1152 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_alltoall_nocomm.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> Sends data from all to all processes, non-mpi version
subroutine mpp_alltoall_int4(sbuf, scount, rbuf, rcount, pelist)
    integer(i4_kind), dimension(:), intent(in) :: sbuf
    integer(i4_kind), dimension(:), intent(inout) :: rbuf
    integer,   intent(in) :: scount, rcount

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call system_clock_default(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, 4)

end subroutine mpp_alltoall_int4


!> Sends data from all to all processes with vector displacement, non-mpi version
subroutine mpp_alltoall_int4_v(sbuf, ssize, sdispl, rbuf, rsize, rdispl, pelist)
    integer(i4_kind), intent(in) :: sbuf(:)
    integer(i4_kind), intent(inout) :: rbuf(:)

    integer, intent(in) :: ssize(:), rsize(:)
    integer, intent(in) :: sdispl(:), rdispl(:)

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call system_clock_default(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, 4)

end subroutine mpp_alltoall_int4_v


!> Sends data from all to all processes with given data types,
!! displacements and block sizes. Non-mpi version
subroutine mpp_alltoall_int4_w(sbuf, ssize, sdispl, stype, &
                          rbuf, rsize, rdispl, rtype, pelist)
    integer(i4_kind), intent(in) :: sbuf(:)
    integer(i4_kind), intent(inout) :: rbuf(:)

    integer, intent(in) :: ssize(:), rsize(:)
    integer, intent(in) :: sdispl(:), rdispl(:)
    type(mpp_type), intent(in) :: stype(:), rtype(:)

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call system_clock_default(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, 4)

end subroutine mpp_alltoall_int4_w
# 1152 "mpp/include/mpp_comm_nocomm.inc" 2

# 1166 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_alltoall_nocomm.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> Sends data from all to all processes, non-mpi version
subroutine mpp_alltoall_int8(sbuf, scount, rbuf, rcount, pelist)
    integer(i8_kind), dimension(:), intent(in) :: sbuf
    integer(i8_kind), dimension(:), intent(inout) :: rbuf
    integer,   intent(in) :: scount, rcount

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call system_clock_default(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, 8)

end subroutine mpp_alltoall_int8


!> Sends data from all to all processes with vector displacement, non-mpi version
subroutine mpp_alltoall_int8_v(sbuf, ssize, sdispl, rbuf, rsize, rdispl, pelist)
    integer(i8_kind), intent(in) :: sbuf(:)
    integer(i8_kind), intent(inout) :: rbuf(:)

    integer, intent(in) :: ssize(:), rsize(:)
    integer, intent(in) :: sdispl(:), rdispl(:)

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call system_clock_default(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, 8)

end subroutine mpp_alltoall_int8_v


!> Sends data from all to all processes with given data types,
!! displacements and block sizes. Non-mpi version
subroutine mpp_alltoall_int8_w(sbuf, ssize, sdispl, stype, &
                          rbuf, rsize, rdispl, rtype, pelist)
    integer(i8_kind), intent(in) :: sbuf(:)
    integer(i8_kind), intent(inout) :: rbuf(:)

    integer, intent(in) :: ssize(:), rsize(:)
    integer, intent(in) :: sdispl(:), rdispl(:)
    type(mpp_type), intent(in) :: stype(:), rtype(:)

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call system_clock_default(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, 8)

end subroutine mpp_alltoall_int8_w
# 1166 "mpp/include/mpp_comm_nocomm.inc" 2

# 1180 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_alltoall_nocomm.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> Sends data from all to all processes, non-mpi version
subroutine mpp_alltoall_real4(sbuf, scount, rbuf, rcount, pelist)
    real(r4_kind), dimension(:), intent(in) :: sbuf
    real(r4_kind), dimension(:), intent(inout) :: rbuf
    integer,   intent(in) :: scount, rcount

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call system_clock_default(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, 4)

end subroutine mpp_alltoall_real4


!> Sends data from all to all processes with vector displacement, non-mpi version
subroutine mpp_alltoall_real4_v(sbuf, ssize, sdispl, rbuf, rsize, rdispl, pelist)
    real(r4_kind), intent(in) :: sbuf(:)
    real(r4_kind), intent(inout) :: rbuf(:)

    integer, intent(in) :: ssize(:), rsize(:)
    integer, intent(in) :: sdispl(:), rdispl(:)

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call system_clock_default(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, 4)

end subroutine mpp_alltoall_real4_v


!> Sends data from all to all processes with given data types,
!! displacements and block sizes. Non-mpi version
subroutine mpp_alltoall_real4_w(sbuf, ssize, sdispl, stype, &
                          rbuf, rsize, rdispl, rtype, pelist)
    real(r4_kind), intent(in) :: sbuf(:)
    real(r4_kind), intent(inout) :: rbuf(:)

    integer, intent(in) :: ssize(:), rsize(:)
    integer, intent(in) :: sdispl(:), rdispl(:)
    type(mpp_type), intent(in) :: stype(:), rtype(:)

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call system_clock_default(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, 4)

end subroutine mpp_alltoall_real4_w
# 1180 "mpp/include/mpp_comm_nocomm.inc" 2

# 1194 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_alltoall_nocomm.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> Sends data from all to all processes, non-mpi version
subroutine mpp_alltoall_real8(sbuf, scount, rbuf, rcount, pelist)
    real(r8_kind), dimension(:), intent(in) :: sbuf
    real(r8_kind), dimension(:), intent(inout) :: rbuf
    integer,   intent(in) :: scount, rcount

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call system_clock_default(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, 8)

end subroutine mpp_alltoall_real8


!> Sends data from all to all processes with vector displacement, non-mpi version
subroutine mpp_alltoall_real8_v(sbuf, ssize, sdispl, rbuf, rsize, rdispl, pelist)
    real(r8_kind), intent(in) :: sbuf(:)
    real(r8_kind), intent(inout) :: rbuf(:)

    integer, intent(in) :: ssize(:), rsize(:)
    integer, intent(in) :: sdispl(:), rdispl(:)

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call system_clock_default(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, 8)

end subroutine mpp_alltoall_real8_v


!> Sends data from all to all processes with given data types,
!! displacements and block sizes. Non-mpi version
subroutine mpp_alltoall_real8_w(sbuf, ssize, sdispl, stype, &
                          rbuf, rsize, rdispl, rtype, pelist)
    real(r8_kind), intent(in) :: sbuf(:)
    real(r8_kind), intent(inout) :: rbuf(:)

    integer, intent(in) :: ssize(:), rsize(:)
    integer, intent(in) :: sdispl(:), rdispl(:)
    type(mpp_type), intent(in) :: stype(:), rtype(:)

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call system_clock_default(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, 8)

end subroutine mpp_alltoall_real8_w
# 1194 "mpp/include/mpp_comm_nocomm.inc" 2

# 1208 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_alltoall_nocomm.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> Sends data from all to all processes, non-mpi version
subroutine mpp_alltoall_logical4(sbuf, scount, rbuf, rcount, pelist)
    logical(l4_kind), dimension(:), intent(in) :: sbuf
    logical(l4_kind), dimension(:), intent(inout) :: rbuf
    integer,   intent(in) :: scount, rcount

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call system_clock_default(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, 4)

end subroutine mpp_alltoall_logical4


!> Sends data from all to all processes with vector displacement, non-mpi version
subroutine mpp_alltoall_logical4_v(sbuf, ssize, sdispl, rbuf, rsize, rdispl, pelist)
    logical(l4_kind), intent(in) :: sbuf(:)
    logical(l4_kind), intent(inout) :: rbuf(:)

    integer, intent(in) :: ssize(:), rsize(:)
    integer, intent(in) :: sdispl(:), rdispl(:)

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call system_clock_default(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, 4)

end subroutine mpp_alltoall_logical4_v


!> Sends data from all to all processes with given data types,
!! displacements and block sizes. Non-mpi version
subroutine mpp_alltoall_logical4_w(sbuf, ssize, sdispl, stype, &
                          rbuf, rsize, rdispl, rtype, pelist)
    logical(l4_kind), intent(in) :: sbuf(:)
    logical(l4_kind), intent(inout) :: rbuf(:)

    integer, intent(in) :: ssize(:), rsize(:)
    integer, intent(in) :: sdispl(:), rdispl(:)
    type(mpp_type), intent(in) :: stype(:), rtype(:)

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call system_clock_default(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, 4)

end subroutine mpp_alltoall_logical4_w
# 1208 "mpp/include/mpp_comm_nocomm.inc" 2

# 1222 "mpp/include/mpp_comm_nocomm.inc"
# 1 "mpp/include/mpp_alltoall_nocomm.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> Sends data from all to all processes, non-mpi version
subroutine mpp_alltoall_logical8(sbuf, scount, rbuf, rcount, pelist)
    logical(l8_kind), dimension(:), intent(in) :: sbuf
    logical(l8_kind), dimension(:), intent(inout) :: rbuf
    integer,   intent(in) :: scount, rcount

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call system_clock_default(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, 8)

end subroutine mpp_alltoall_logical8


!> Sends data from all to all processes with vector displacement, non-mpi version
subroutine mpp_alltoall_logical8_v(sbuf, ssize, sdispl, rbuf, rsize, rdispl, pelist)
    logical(l8_kind), intent(in) :: sbuf(:)
    logical(l8_kind), intent(inout) :: rbuf(:)

    integer, intent(in) :: ssize(:), rsize(:)
    integer, intent(in) :: sdispl(:), rdispl(:)

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call system_clock_default(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, 8)

end subroutine mpp_alltoall_logical8_v


!> Sends data from all to all processes with given data types,
!! displacements and block sizes. Non-mpi version
subroutine mpp_alltoall_logical8_w(sbuf, ssize, sdispl, stype, &
                          rbuf, rsize, rdispl, rtype, pelist)
    logical(l8_kind), intent(in) :: sbuf(:)
    logical(l8_kind), intent(inout) :: rbuf(:)

    integer, intent(in) :: ssize(:), rsize(:)
    integer, intent(in) :: sdispl(:), rdispl(:)
    type(mpp_type), intent(in) :: stype(:), rtype(:)

    integer, intent(in), optional :: pelist(0:)

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_ALLTOALL: You must first call mpp_init.')

    if (current_clock .NE. 0) call system_clock_default(start_tick)

    rbuf(:) = sbuf(:)

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_ALLTOALL, 8)

end subroutine mpp_alltoall_logical8_w
# 1222 "mpp/include/mpp_comm_nocomm.inc" 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!            DATA TRANSFER TYPES: mpp_type_create, mpp_type_free              !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





# 1 "mpp/include/mpp_type_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

subroutine mpp_type_create_int4(field, array_of_subsizes, array_of_starts, dtype)
    integer(i4_kind), intent(in) :: field(:,:,:)
    integer, intent(in) :: array_of_subsizes(:)
    integer, intent(in) :: array_of_starts(:)
    type(mpp_type), target, intent(out) :: dtype

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_TYPE_CREATE: You must first call mpp_init.')

    if (current_clock .NE. 0) &
        call system_clock_default(start_tick)

    call mpp_error(NOTE, 'MPP_TYPE_CREATE: &
                         &This function is not used in serial mode.')

    ! For consistency with the MPI version, we return a valid mpp_type
    dtype = mpp_byte

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_TYPE_CREATE, 8)

end subroutine mpp_type_create_int4
# 1232 "mpp/include/mpp_comm_nocomm.inc" 2





# 1 "mpp/include/mpp_type_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

subroutine mpp_type_create_int8(field, array_of_subsizes, array_of_starts, dtype)
    integer(i8_kind), intent(in) :: field(:,:,:)
    integer, intent(in) :: array_of_subsizes(:)
    integer, intent(in) :: array_of_starts(:)
    type(mpp_type), target, intent(out) :: dtype

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_TYPE_CREATE: You must first call mpp_init.')

    if (current_clock .NE. 0) &
        call system_clock_default(start_tick)

    call mpp_error(NOTE, 'MPP_TYPE_CREATE: &
                         &This function is not used in serial mode.')

    ! For consistency with the MPI version, we return a valid mpp_type
    dtype = mpp_byte

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_TYPE_CREATE, 8)

end subroutine mpp_type_create_int8
# 1237 "mpp/include/mpp_comm_nocomm.inc" 2





# 1 "mpp/include/mpp_type_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

subroutine mpp_type_create_real4(field, array_of_subsizes, array_of_starts, dtype)
    real(r4_kind), intent(in) :: field(:,:,:)
    integer, intent(in) :: array_of_subsizes(:)
    integer, intent(in) :: array_of_starts(:)
    type(mpp_type), target, intent(out) :: dtype

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_TYPE_CREATE: You must first call mpp_init.')

    if (current_clock .NE. 0) &
        call system_clock_default(start_tick)

    call mpp_error(NOTE, 'MPP_TYPE_CREATE: &
                         &This function is not used in serial mode.')

    ! For consistency with the MPI version, we return a valid mpp_type
    dtype = mpp_byte

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_TYPE_CREATE, 8)

end subroutine mpp_type_create_real4
# 1242 "mpp/include/mpp_comm_nocomm.inc" 2





# 1 "mpp/include/mpp_type_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

subroutine mpp_type_create_real8(field, array_of_subsizes, array_of_starts, dtype)
    real(r8_kind), intent(in) :: field(:,:,:)
    integer, intent(in) :: array_of_subsizes(:)
    integer, intent(in) :: array_of_starts(:)
    type(mpp_type), target, intent(out) :: dtype

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_TYPE_CREATE: You must first call mpp_init.')

    if (current_clock .NE. 0) &
        call system_clock_default(start_tick)

    call mpp_error(NOTE, 'MPP_TYPE_CREATE: &
                         &This function is not used in serial mode.')

    ! For consistency with the MPI version, we return a valid mpp_type
    dtype = mpp_byte

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_TYPE_CREATE, 8)

end subroutine mpp_type_create_real8
# 1247 "mpp/include/mpp_comm_nocomm.inc" 2





# 1 "mpp/include/mpp_type_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

subroutine mpp_type_create_cmplx4(field, array_of_subsizes, array_of_starts, dtype)
    complex(c4_kind), intent(in) :: field(:,:,:)
    integer, intent(in) :: array_of_subsizes(:)
    integer, intent(in) :: array_of_starts(:)
    type(mpp_type), target, intent(out) :: dtype

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_TYPE_CREATE: You must first call mpp_init.')

    if (current_clock .NE. 0) &
        call system_clock_default(start_tick)

    call mpp_error(NOTE, 'MPP_TYPE_CREATE: &
                         &This function is not used in serial mode.')

    ! For consistency with the MPI version, we return a valid mpp_type
    dtype = mpp_byte

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_TYPE_CREATE, 8)

end subroutine mpp_type_create_cmplx4
# 1252 "mpp/include/mpp_comm_nocomm.inc" 2





# 1 "mpp/include/mpp_type_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

subroutine mpp_type_create_cmplx8(field, array_of_subsizes, array_of_starts, dtype)
    complex(c8_kind), intent(in) :: field(:,:,:)
    integer, intent(in) :: array_of_subsizes(:)
    integer, intent(in) :: array_of_starts(:)
    type(mpp_type), target, intent(out) :: dtype

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_TYPE_CREATE: You must first call mpp_init.')

    if (current_clock .NE. 0) &
        call system_clock_default(start_tick)

    call mpp_error(NOTE, 'MPP_TYPE_CREATE: &
                         &This function is not used in serial mode.')

    ! For consistency with the MPI version, we return a valid mpp_type
    dtype = mpp_byte

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_TYPE_CREATE, 8)

end subroutine mpp_type_create_cmplx8
# 1257 "mpp/include/mpp_comm_nocomm.inc" 2





# 1 "mpp/include/mpp_type_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

subroutine mpp_type_create_logical4(field, array_of_subsizes, array_of_starts, dtype)
    logical(l4_kind), intent(in) :: field(:,:,:)
    integer, intent(in) :: array_of_subsizes(:)
    integer, intent(in) :: array_of_starts(:)
    type(mpp_type), target, intent(out) :: dtype

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_TYPE_CREATE: You must first call mpp_init.')

    if (current_clock .NE. 0) &
        call system_clock_default(start_tick)

    call mpp_error(NOTE, 'MPP_TYPE_CREATE: &
                         &This function is not used in serial mode.')

    ! For consistency with the MPI version, we return a valid mpp_type
    dtype = mpp_byte

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_TYPE_CREATE, 8)

end subroutine mpp_type_create_logical4
# 1262 "mpp/include/mpp_comm_nocomm.inc" 2





# 1 "mpp/include/mpp_type_nocomm.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

subroutine mpp_type_create_logical8(field, array_of_subsizes, array_of_starts, dtype)
    logical(l8_kind), intent(in) :: field(:,:,:)
    integer, intent(in) :: array_of_subsizes(:)
    integer, intent(in) :: array_of_starts(:)
    type(mpp_type), target, intent(out) :: dtype

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_TYPE_CREATE: You must first call mpp_init.')

    if (current_clock .NE. 0) &
        call system_clock_default(start_tick)

    call mpp_error(NOTE, 'MPP_TYPE_CREATE: &
                         &This function is not used in serial mode.')

    ! For consistency with the MPI version, we return a valid mpp_type
    dtype = mpp_byte

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_TYPE_CREATE, 8)

end subroutine mpp_type_create_logical8
# 1267 "mpp/include/mpp_comm_nocomm.inc" 2

! Clear preprocessor flags




subroutine mpp_type_free(dtype)
    type(mpp_type), pointer, intent(inout) :: dtype

    call mpp_error(NOTE, 'MPP_TYPE_FREE: ' &
                         //'This function should not be used in serial mode.')

    ! For consistency with MPI, we deallocate the pointer
    deallocate(dtype)

end subroutine mpp_type_free
# 38 "mpp/include/mpp_comm.inc" 2


# 49 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum_int.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> Calculates integer checksum over pelist
function mpp_chksum_i8_1d( var, pelist, mask_val )
  integer(i8_kind) :: mpp_chksum_i8_1d
  integer(i8_kind), intent(in) :: var (:)
  integer, optional :: pelist(:)
  integer(i8_kind), intent(in), optional :: mask_val

  if ( PRESENT(mask_val) ) then
     !PACK on var/=mask_val ignores values in var
     !equiv to setting those values=0, but on sparse arrays
     !pack should return much smaller array to sum
     mpp_chksum_i8_1d = sum( INT( PACK(var,var/=mask_val),i8_kind) )
  else
     mpp_chksum_i8_1d = sum(INT(var,i8_kind))
  end if

      call mpp_sum( mpp_chksum_i8_1d, pelist )
      return

    end function mpp_chksum_i8_1d


!> Handles real mask for easier implementation
function mpp_chksum_i8_1d_rmask( var, pelist, mask_val )
  integer(KIND=i8_kind) :: mpp_chksum_i8_1d_rmask
  integer(i8_kind), intent(in) :: var (:)
  integer, optional :: pelist(:)
  real, intent(in) :: mask_val
  integer(KIND(var))::imask_val
  integer(KIND=i4_kind)::i4tmp(2)=0
  real(KIND=r4_kind)::r4tmp(2)=0
  integer(KIND=i8_kind) :: i8tmp=0
  !high fidelity error message
  character(LEN=1) :: tmpStr1,tmpStr2,tmpStr3
  character(LEN=32) :: tmpStr4
  character(LEN=512) :: errStr

! Primary Logic: These first two are the "expected" branches.
!! These all resolve to MPP_FILL_INT
  !!Should catch real "default_fill"(MPP_FILL_DOUBLE)
  if (mask_val == MPP_FILL_DOUBLE ) then !this is FMS variable field default fill
     ! we've packed an MPP_FILL_
     imask_val = MPP_FILL_INT
 !Current NETCDF fill values (AKA MPP_FILL_*) designed towards CEILING(MPP_FILL_{FLOAT,DOUBLE},kind=4byte)=MPP_FILL_INT
  else if ( CEILING(mask_val, i4_kind) == MPP_FILL_INT ) then
     ! we've also packed an MPP_FILL_
     imask_val = MPP_FILL_INT
! Secondary Logic:
!! We've done something dangerous
  else
     i8tmp = TRANSFER(mask_val , i8tmp )
     i4tmp = TRANSFER(mask_val , i4tmp )
     r4tmp = TRANSFER(mask_val , r4tmp )
     if ( i8tmp == MPP_FILL_INT ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(i4tmp == MPP_FILL_INT) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(r4tmp == MPP_FILL_DOUBLE) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else
        ! we have no idea what this is
        ! construct detailed errStr
        errStr = "mpp_chksum: mpp_chksum_i"
        write(unit=tmpStr1,fmt="(I1)") KIND(var)
        write(unit=tmpstr2,fmt="(I1)") SIZE(SHAPE(var))
        errStr = errStr // tmpStr1 // "_" // tmpstr2 // "d_rmask passed int var with REAL("
        write(unit=tmpstr3,fmt="(I1)") KIND(mask_val)
        errStr = errStr // tmpstr3 // ") mask_val="
        write(unit=tmpstr4,fmt=*) mask_val
        errStr = errStr // trim(tmpstr4) // "has been called with these strange values. Check your KINDS, "// &
              "_FillValue, pack and mask_val. "// &
              "Hint: Try being explicit and using MPP_FILL_{INT,FLOAT,DOUBLE}. "// &
              "Continuing by using the default MPP_FILL_INT. " // &
              "THIS WILL BE FATAL IN THE FUTURE!"
        call mpp_error(WARNING, trim(errStr) )

        imask_val = MPP_FILL_INT
     end if
  end if

  mpp_chksum_i8_1d_rmask = mpp_chksum(var,pelist,mask_val=imask_val)

  return

end function mpp_chksum_i8_1d_rmask
# 49 "mpp/include/mpp_comm.inc" 2

# 59 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum_int.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> Calculates integer checksum over pelist
function mpp_chksum_i8_2d( var, pelist, mask_val )
  integer(i8_kind) :: mpp_chksum_i8_2d
  integer(i8_kind), intent(in) :: var (:,:)
  integer, optional :: pelist(:)
  integer(i8_kind), intent(in), optional :: mask_val

  if ( PRESENT(mask_val) ) then
     !PACK on var/=mask_val ignores values in var
     !equiv to setting those values=0, but on sparse arrays
     !pack should return much smaller array to sum
     mpp_chksum_i8_2d = sum( INT( PACK(var,var/=mask_val),i8_kind) )
  else
     mpp_chksum_i8_2d = sum(INT(var,i8_kind))
  end if

      call mpp_sum( mpp_chksum_i8_2d, pelist )
      return

    end function mpp_chksum_i8_2d


!> Handles real mask for easier implementation
function mpp_chksum_i8_2d_rmask( var, pelist, mask_val )
  integer(KIND=i8_kind) :: mpp_chksum_i8_2d_rmask
  integer(i8_kind), intent(in) :: var (:,:)
  integer, optional :: pelist(:)
  real, intent(in) :: mask_val
  integer(KIND(var))::imask_val
  integer(KIND=i4_kind)::i4tmp(2)=0
  real(KIND=r4_kind)::r4tmp(2)=0
  integer(KIND=i8_kind) :: i8tmp=0
  !high fidelity error message
  character(LEN=1) :: tmpStr1,tmpStr2,tmpStr3
  character(LEN=32) :: tmpStr4
  character(LEN=512) :: errStr

! Primary Logic: These first two are the "expected" branches.
!! These all resolve to MPP_FILL_INT
  !!Should catch real "default_fill"(MPP_FILL_DOUBLE)
  if (mask_val == MPP_FILL_DOUBLE ) then !this is FMS variable field default fill
     ! we've packed an MPP_FILL_
     imask_val = MPP_FILL_INT
 !Current NETCDF fill values (AKA MPP_FILL_*) designed towards CEILING(MPP_FILL_{FLOAT,DOUBLE},kind=4byte)=MPP_FILL_INT
  else if ( CEILING(mask_val, i4_kind) == MPP_FILL_INT ) then
     ! we've also packed an MPP_FILL_
     imask_val = MPP_FILL_INT
! Secondary Logic:
!! We've done something dangerous
  else
     i8tmp = TRANSFER(mask_val , i8tmp )
     i4tmp = TRANSFER(mask_val , i4tmp )
     r4tmp = TRANSFER(mask_val , r4tmp )
     if ( i8tmp == MPP_FILL_INT ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(i4tmp == MPP_FILL_INT) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(r4tmp == MPP_FILL_DOUBLE) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else
        ! we have no idea what this is
        ! construct detailed errStr
        errStr = "mpp_chksum: mpp_chksum_i"
        write(unit=tmpStr1,fmt="(I1)") KIND(var)
        write(unit=tmpstr2,fmt="(I1)") SIZE(SHAPE(var))
        errStr = errStr // tmpStr1 // "_" // tmpstr2 // "d_rmask passed int var with REAL("
        write(unit=tmpstr3,fmt="(I1)") KIND(mask_val)
        errStr = errStr // tmpstr3 // ") mask_val="
        write(unit=tmpstr4,fmt=*) mask_val
        errStr = errStr // trim(tmpstr4) // "has been called with these strange values. Check your KINDS, "// &
              "_FillValue, pack and mask_val. "// &
              "Hint: Try being explicit and using MPP_FILL_{INT,FLOAT,DOUBLE}. "// &
              "Continuing by using the default MPP_FILL_INT. " // &
              "THIS WILL BE FATAL IN THE FUTURE!"
        call mpp_error(WARNING, trim(errStr) )

        imask_val = MPP_FILL_INT
     end if
  end if

  mpp_chksum_i8_2d_rmask = mpp_chksum(var,pelist,mask_val=imask_val)

  return

end function mpp_chksum_i8_2d_rmask
# 59 "mpp/include/mpp_comm.inc" 2

# 69 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum_int.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> Calculates integer checksum over pelist
function mpp_chksum_i8_3d( var, pelist, mask_val )
  integer(i8_kind) :: mpp_chksum_i8_3d
  integer(i8_kind), intent(in) :: var (:,:,:)
  integer, optional :: pelist(:)
  integer(i8_kind), intent(in), optional :: mask_val

  if ( PRESENT(mask_val) ) then
     !PACK on var/=mask_val ignores values in var
     !equiv to setting those values=0, but on sparse arrays
     !pack should return much smaller array to sum
     mpp_chksum_i8_3d = sum( INT( PACK(var,var/=mask_val),i8_kind) )
  else
     mpp_chksum_i8_3d = sum(INT(var,i8_kind))
  end if

      call mpp_sum( mpp_chksum_i8_3d, pelist )
      return

    end function mpp_chksum_i8_3d


!> Handles real mask for easier implementation
function mpp_chksum_i8_3d_rmask( var, pelist, mask_val )
  integer(KIND=i8_kind) :: mpp_chksum_i8_3d_rmask
  integer(i8_kind), intent(in) :: var (:,:,:)
  integer, optional :: pelist(:)
  real, intent(in) :: mask_val
  integer(KIND(var))::imask_val
  integer(KIND=i4_kind)::i4tmp(2)=0
  real(KIND=r4_kind)::r4tmp(2)=0
  integer(KIND=i8_kind) :: i8tmp=0
  !high fidelity error message
  character(LEN=1) :: tmpStr1,tmpStr2,tmpStr3
  character(LEN=32) :: tmpStr4
  character(LEN=512) :: errStr

! Primary Logic: These first two are the "expected" branches.
!! These all resolve to MPP_FILL_INT
  !!Should catch real "default_fill"(MPP_FILL_DOUBLE)
  if (mask_val == MPP_FILL_DOUBLE ) then !this is FMS variable field default fill
     ! we've packed an MPP_FILL_
     imask_val = MPP_FILL_INT
 !Current NETCDF fill values (AKA MPP_FILL_*) designed towards CEILING(MPP_FILL_{FLOAT,DOUBLE},kind=4byte)=MPP_FILL_INT
  else if ( CEILING(mask_val, i4_kind) == MPP_FILL_INT ) then
     ! we've also packed an MPP_FILL_
     imask_val = MPP_FILL_INT
! Secondary Logic:
!! We've done something dangerous
  else
     i8tmp = TRANSFER(mask_val , i8tmp )
     i4tmp = TRANSFER(mask_val , i4tmp )
     r4tmp = TRANSFER(mask_val , r4tmp )
     if ( i8tmp == MPP_FILL_INT ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(i4tmp == MPP_FILL_INT) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(r4tmp == MPP_FILL_DOUBLE) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else
        ! we have no idea what this is
        ! construct detailed errStr
        errStr = "mpp_chksum: mpp_chksum_i"
        write(unit=tmpStr1,fmt="(I1)") KIND(var)
        write(unit=tmpstr2,fmt="(I1)") SIZE(SHAPE(var))
        errStr = errStr // tmpStr1 // "_" // tmpstr2 // "d_rmask passed int var with REAL("
        write(unit=tmpstr3,fmt="(I1)") KIND(mask_val)
        errStr = errStr // tmpstr3 // ") mask_val="
        write(unit=tmpstr4,fmt=*) mask_val
        errStr = errStr // trim(tmpstr4) // "has been called with these strange values. Check your KINDS, "// &
              "_FillValue, pack and mask_val. "// &
              "Hint: Try being explicit and using MPP_FILL_{INT,FLOAT,DOUBLE}. "// &
              "Continuing by using the default MPP_FILL_INT. " // &
              "THIS WILL BE FATAL IN THE FUTURE!"
        call mpp_error(WARNING, trim(errStr) )

        imask_val = MPP_FILL_INT
     end if
  end if

  mpp_chksum_i8_3d_rmask = mpp_chksum(var,pelist,mask_val=imask_val)

  return

end function mpp_chksum_i8_3d_rmask
# 69 "mpp/include/mpp_comm.inc" 2

# 79 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum_int.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> Calculates integer checksum over pelist
function mpp_chksum_i8_4d( var, pelist, mask_val )
  integer(i8_kind) :: mpp_chksum_i8_4d
  integer(i8_kind), intent(in) :: var (:,:,:,:)
  integer, optional :: pelist(:)
  integer(i8_kind), intent(in), optional :: mask_val

  if ( PRESENT(mask_val) ) then
     !PACK on var/=mask_val ignores values in var
     !equiv to setting those values=0, but on sparse arrays
     !pack should return much smaller array to sum
     mpp_chksum_i8_4d = sum( INT( PACK(var,var/=mask_val),i8_kind) )
  else
     mpp_chksum_i8_4d = sum(INT(var,i8_kind))
  end if

      call mpp_sum( mpp_chksum_i8_4d, pelist )
      return

    end function mpp_chksum_i8_4d


!> Handles real mask for easier implementation
function mpp_chksum_i8_4d_rmask( var, pelist, mask_val )
  integer(KIND=i8_kind) :: mpp_chksum_i8_4d_rmask
  integer(i8_kind), intent(in) :: var (:,:,:,:)
  integer, optional :: pelist(:)
  real, intent(in) :: mask_val
  integer(KIND(var))::imask_val
  integer(KIND=i4_kind)::i4tmp(2)=0
  real(KIND=r4_kind)::r4tmp(2)=0
  integer(KIND=i8_kind) :: i8tmp=0
  !high fidelity error message
  character(LEN=1) :: tmpStr1,tmpStr2,tmpStr3
  character(LEN=32) :: tmpStr4
  character(LEN=512) :: errStr

! Primary Logic: These first two are the "expected" branches.
!! These all resolve to MPP_FILL_INT
  !!Should catch real "default_fill"(MPP_FILL_DOUBLE)
  if (mask_val == MPP_FILL_DOUBLE ) then !this is FMS variable field default fill
     ! we've packed an MPP_FILL_
     imask_val = MPP_FILL_INT
 !Current NETCDF fill values (AKA MPP_FILL_*) designed towards CEILING(MPP_FILL_{FLOAT,DOUBLE},kind=4byte)=MPP_FILL_INT
  else if ( CEILING(mask_val, i4_kind) == MPP_FILL_INT ) then
     ! we've also packed an MPP_FILL_
     imask_val = MPP_FILL_INT
! Secondary Logic:
!! We've done something dangerous
  else
     i8tmp = TRANSFER(mask_val , i8tmp )
     i4tmp = TRANSFER(mask_val , i4tmp )
     r4tmp = TRANSFER(mask_val , r4tmp )
     if ( i8tmp == MPP_FILL_INT ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(i4tmp == MPP_FILL_INT) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(r4tmp == MPP_FILL_DOUBLE) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else
        ! we have no idea what this is
        ! construct detailed errStr
        errStr = "mpp_chksum: mpp_chksum_i"
        write(unit=tmpStr1,fmt="(I1)") KIND(var)
        write(unit=tmpstr2,fmt="(I1)") SIZE(SHAPE(var))
        errStr = errStr // tmpStr1 // "_" // tmpstr2 // "d_rmask passed int var with REAL("
        write(unit=tmpstr3,fmt="(I1)") KIND(mask_val)
        errStr = errStr // tmpstr3 // ") mask_val="
        write(unit=tmpstr4,fmt=*) mask_val
        errStr = errStr // trim(tmpstr4) // "has been called with these strange values. Check your KINDS, "// &
              "_FillValue, pack and mask_val. "// &
              "Hint: Try being explicit and using MPP_FILL_{INT,FLOAT,DOUBLE}. "// &
              "Continuing by using the default MPP_FILL_INT. " // &
              "THIS WILL BE FATAL IN THE FUTURE!"
        call mpp_error(WARNING, trim(errStr) )

        imask_val = MPP_FILL_INT
     end if
  end if

  mpp_chksum_i8_4d_rmask = mpp_chksum(var,pelist,mask_val=imask_val)

  return

end function mpp_chksum_i8_4d_rmask
# 79 "mpp/include/mpp_comm.inc" 2

# 89 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum_int.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> Calculates integer checksum over pelist
function mpp_chksum_i8_5d( var, pelist, mask_val )
  integer(i8_kind) :: mpp_chksum_i8_5d
  integer(i8_kind), intent(in) :: var (:,:,:,:,:)
  integer, optional :: pelist(:)
  integer(i8_kind), intent(in), optional :: mask_val

  if ( PRESENT(mask_val) ) then
     !PACK on var/=mask_val ignores values in var
     !equiv to setting those values=0, but on sparse arrays
     !pack should return much smaller array to sum
     mpp_chksum_i8_5d = sum( INT( PACK(var,var/=mask_val),i8_kind) )
  else
     mpp_chksum_i8_5d = sum(INT(var,i8_kind))
  end if

      call mpp_sum( mpp_chksum_i8_5d, pelist )
      return

    end function mpp_chksum_i8_5d


!> Handles real mask for easier implementation
function mpp_chksum_i8_5d_rmask( var, pelist, mask_val )
  integer(KIND=i8_kind) :: mpp_chksum_i8_5d_rmask
  integer(i8_kind), intent(in) :: var (:,:,:,:,:)
  integer, optional :: pelist(:)
  real, intent(in) :: mask_val
  integer(KIND(var))::imask_val
  integer(KIND=i4_kind)::i4tmp(2)=0
  real(KIND=r4_kind)::r4tmp(2)=0
  integer(KIND=i8_kind) :: i8tmp=0
  !high fidelity error message
  character(LEN=1) :: tmpStr1,tmpStr2,tmpStr3
  character(LEN=32) :: tmpStr4
  character(LEN=512) :: errStr

! Primary Logic: These first two are the "expected" branches.
!! These all resolve to MPP_FILL_INT
  !!Should catch real "default_fill"(MPP_FILL_DOUBLE)
  if (mask_val == MPP_FILL_DOUBLE ) then !this is FMS variable field default fill
     ! we've packed an MPP_FILL_
     imask_val = MPP_FILL_INT
 !Current NETCDF fill values (AKA MPP_FILL_*) designed towards CEILING(MPP_FILL_{FLOAT,DOUBLE},kind=4byte)=MPP_FILL_INT
  else if ( CEILING(mask_val, i4_kind) == MPP_FILL_INT ) then
     ! we've also packed an MPP_FILL_
     imask_val = MPP_FILL_INT
! Secondary Logic:
!! We've done something dangerous
  else
     i8tmp = TRANSFER(mask_val , i8tmp )
     i4tmp = TRANSFER(mask_val , i4tmp )
     r4tmp = TRANSFER(mask_val , r4tmp )
     if ( i8tmp == MPP_FILL_INT ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(i4tmp == MPP_FILL_INT) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(r4tmp == MPP_FILL_DOUBLE) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else
        ! we have no idea what this is
        ! construct detailed errStr
        errStr = "mpp_chksum: mpp_chksum_i"
        write(unit=tmpStr1,fmt="(I1)") KIND(var)
        write(unit=tmpstr2,fmt="(I1)") SIZE(SHAPE(var))
        errStr = errStr // tmpStr1 // "_" // tmpstr2 // "d_rmask passed int var with REAL("
        write(unit=tmpstr3,fmt="(I1)") KIND(mask_val)
        errStr = errStr // tmpstr3 // ") mask_val="
        write(unit=tmpstr4,fmt=*) mask_val
        errStr = errStr // trim(tmpstr4) // "has been called with these strange values. Check your KINDS, "// &
              "_FillValue, pack and mask_val. "// &
              "Hint: Try being explicit and using MPP_FILL_{INT,FLOAT,DOUBLE}. "// &
              "Continuing by using the default MPP_FILL_INT. " // &
              "THIS WILL BE FATAL IN THE FUTURE!"
        call mpp_error(WARNING, trim(errStr) )

        imask_val = MPP_FILL_INT
     end if
  end if

  mpp_chksum_i8_5d_rmask = mpp_chksum(var,pelist,mask_val=imask_val)

  return

end function mpp_chksum_i8_5d_rmask
# 89 "mpp/include/mpp_comm.inc" 2

# 99 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum_int.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> Calculates integer checksum over pelist
function mpp_chksum_i4_1d( var, pelist, mask_val )
  integer(i8_kind) :: mpp_chksum_i4_1d
  integer(i4_kind), intent(in) :: var (:)
  integer, optional :: pelist(:)
  integer(i4_kind), intent(in), optional :: mask_val

  if ( PRESENT(mask_val) ) then
     !PACK on var/=mask_val ignores values in var
     !equiv to setting those values=0, but on sparse arrays
     !pack should return much smaller array to sum
     mpp_chksum_i4_1d = sum( INT( PACK(var,var/=mask_val),i8_kind) )
  else
     mpp_chksum_i4_1d = sum(INT(var,i8_kind))
  end if

      call mpp_sum( mpp_chksum_i4_1d, pelist )
      return

    end function mpp_chksum_i4_1d


!> Handles real mask for easier implementation
function mpp_chksum_i4_1d_rmask( var, pelist, mask_val )
  integer(KIND=i8_kind) :: mpp_chksum_i4_1d_rmask
  integer(i4_kind), intent(in) :: var (:)
  integer, optional :: pelist(:)
  real, intent(in) :: mask_val
  integer(KIND(var))::imask_val
  integer(KIND=i4_kind)::i4tmp(2)=0
  real(KIND=r4_kind)::r4tmp(2)=0
  integer(KIND=i8_kind) :: i8tmp=0
  !high fidelity error message
  character(LEN=1) :: tmpStr1,tmpStr2,tmpStr3
  character(LEN=32) :: tmpStr4
  character(LEN=512) :: errStr

! Primary Logic: These first two are the "expected" branches.
!! These all resolve to MPP_FILL_INT
  !!Should catch real "default_fill"(MPP_FILL_DOUBLE)
  if (mask_val == MPP_FILL_DOUBLE ) then !this is FMS variable field default fill
     ! we've packed an MPP_FILL_
     imask_val = MPP_FILL_INT
 !Current NETCDF fill values (AKA MPP_FILL_*) designed towards CEILING(MPP_FILL_{FLOAT,DOUBLE},kind=4byte)=MPP_FILL_INT
  else if ( CEILING(mask_val, i4_kind) == MPP_FILL_INT ) then
     ! we've also packed an MPP_FILL_
     imask_val = MPP_FILL_INT
! Secondary Logic:
!! We've done something dangerous
  else
     i8tmp = TRANSFER(mask_val , i8tmp )
     i4tmp = TRANSFER(mask_val , i4tmp )
     r4tmp = TRANSFER(mask_val , r4tmp )
     if ( i8tmp == MPP_FILL_INT ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(i4tmp == MPP_FILL_INT) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(r4tmp == MPP_FILL_DOUBLE) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else
        ! we have no idea what this is
        ! construct detailed errStr
        errStr = "mpp_chksum: mpp_chksum_i"
        write(unit=tmpStr1,fmt="(I1)") KIND(var)
        write(unit=tmpstr2,fmt="(I1)") SIZE(SHAPE(var))
        errStr = errStr // tmpStr1 // "_" // tmpstr2 // "d_rmask passed int var with REAL("
        write(unit=tmpstr3,fmt="(I1)") KIND(mask_val)
        errStr = errStr // tmpstr3 // ") mask_val="
        write(unit=tmpstr4,fmt=*) mask_val
        errStr = errStr // trim(tmpstr4) // "has been called with these strange values. Check your KINDS, "// &
              "_FillValue, pack and mask_val. "// &
              "Hint: Try being explicit and using MPP_FILL_{INT,FLOAT,DOUBLE}. "// &
              "Continuing by using the default MPP_FILL_INT. " // &
              "THIS WILL BE FATAL IN THE FUTURE!"
        call mpp_error(WARNING, trim(errStr) )

        imask_val = MPP_FILL_INT
     end if
  end if

  mpp_chksum_i4_1d_rmask = mpp_chksum(var,pelist,mask_val=imask_val)

  return

end function mpp_chksum_i4_1d_rmask
# 99 "mpp/include/mpp_comm.inc" 2

# 109 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum_int.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> Calculates integer checksum over pelist
function mpp_chksum_i4_2d( var, pelist, mask_val )
  integer(i8_kind) :: mpp_chksum_i4_2d
  integer(i4_kind), intent(in) :: var (:,:)
  integer, optional :: pelist(:)
  integer(i4_kind), intent(in), optional :: mask_val

  if ( PRESENT(mask_val) ) then
     !PACK on var/=mask_val ignores values in var
     !equiv to setting those values=0, but on sparse arrays
     !pack should return much smaller array to sum
     mpp_chksum_i4_2d = sum( INT( PACK(var,var/=mask_val),i8_kind) )
  else
     mpp_chksum_i4_2d = sum(INT(var,i8_kind))
  end if

      call mpp_sum( mpp_chksum_i4_2d, pelist )
      return

    end function mpp_chksum_i4_2d


!> Handles real mask for easier implementation
function mpp_chksum_i4_2d_rmask( var, pelist, mask_val )
  integer(KIND=i8_kind) :: mpp_chksum_i4_2d_rmask
  integer(i4_kind), intent(in) :: var (:,:)
  integer, optional :: pelist(:)
  real, intent(in) :: mask_val
  integer(KIND(var))::imask_val
  integer(KIND=i4_kind)::i4tmp(2)=0
  real(KIND=r4_kind)::r4tmp(2)=0
  integer(KIND=i8_kind) :: i8tmp=0
  !high fidelity error message
  character(LEN=1) :: tmpStr1,tmpStr2,tmpStr3
  character(LEN=32) :: tmpStr4
  character(LEN=512) :: errStr

! Primary Logic: These first two are the "expected" branches.
!! These all resolve to MPP_FILL_INT
  !!Should catch real "default_fill"(MPP_FILL_DOUBLE)
  if (mask_val == MPP_FILL_DOUBLE ) then !this is FMS variable field default fill
     ! we've packed an MPP_FILL_
     imask_val = MPP_FILL_INT
 !Current NETCDF fill values (AKA MPP_FILL_*) designed towards CEILING(MPP_FILL_{FLOAT,DOUBLE},kind=4byte)=MPP_FILL_INT
  else if ( CEILING(mask_val, i4_kind) == MPP_FILL_INT ) then
     ! we've also packed an MPP_FILL_
     imask_val = MPP_FILL_INT
! Secondary Logic:
!! We've done something dangerous
  else
     i8tmp = TRANSFER(mask_val , i8tmp )
     i4tmp = TRANSFER(mask_val , i4tmp )
     r4tmp = TRANSFER(mask_val , r4tmp )
     if ( i8tmp == MPP_FILL_INT ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(i4tmp == MPP_FILL_INT) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(r4tmp == MPP_FILL_DOUBLE) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else
        ! we have no idea what this is
        ! construct detailed errStr
        errStr = "mpp_chksum: mpp_chksum_i"
        write(unit=tmpStr1,fmt="(I1)") KIND(var)
        write(unit=tmpstr2,fmt="(I1)") SIZE(SHAPE(var))
        errStr = errStr // tmpStr1 // "_" // tmpstr2 // "d_rmask passed int var with REAL("
        write(unit=tmpstr3,fmt="(I1)") KIND(mask_val)
        errStr = errStr // tmpstr3 // ") mask_val="
        write(unit=tmpstr4,fmt=*) mask_val
        errStr = errStr // trim(tmpstr4) // "has been called with these strange values. Check your KINDS, "// &
              "_FillValue, pack and mask_val. "// &
              "Hint: Try being explicit and using MPP_FILL_{INT,FLOAT,DOUBLE}. "// &
              "Continuing by using the default MPP_FILL_INT. " // &
              "THIS WILL BE FATAL IN THE FUTURE!"
        call mpp_error(WARNING, trim(errStr) )

        imask_val = MPP_FILL_INT
     end if
  end if

  mpp_chksum_i4_2d_rmask = mpp_chksum(var,pelist,mask_val=imask_val)

  return

end function mpp_chksum_i4_2d_rmask
# 109 "mpp/include/mpp_comm.inc" 2

# 119 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum_int.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> Calculates integer checksum over pelist
function mpp_chksum_i4_3d( var, pelist, mask_val )
  integer(i8_kind) :: mpp_chksum_i4_3d
  integer(i4_kind), intent(in) :: var (:,:,:)
  integer, optional :: pelist(:)
  integer(i4_kind), intent(in), optional :: mask_val

  if ( PRESENT(mask_val) ) then
     !PACK on var/=mask_val ignores values in var
     !equiv to setting those values=0, but on sparse arrays
     !pack should return much smaller array to sum
     mpp_chksum_i4_3d = sum( INT( PACK(var,var/=mask_val),i8_kind) )
  else
     mpp_chksum_i4_3d = sum(INT(var,i8_kind))
  end if

      call mpp_sum( mpp_chksum_i4_3d, pelist )
      return

    end function mpp_chksum_i4_3d


!> Handles real mask for easier implementation
function mpp_chksum_i4_3d_rmask( var, pelist, mask_val )
  integer(KIND=i8_kind) :: mpp_chksum_i4_3d_rmask
  integer(i4_kind), intent(in) :: var (:,:,:)
  integer, optional :: pelist(:)
  real, intent(in) :: mask_val
  integer(KIND(var))::imask_val
  integer(KIND=i4_kind)::i4tmp(2)=0
  real(KIND=r4_kind)::r4tmp(2)=0
  integer(KIND=i8_kind) :: i8tmp=0
  !high fidelity error message
  character(LEN=1) :: tmpStr1,tmpStr2,tmpStr3
  character(LEN=32) :: tmpStr4
  character(LEN=512) :: errStr

! Primary Logic: These first two are the "expected" branches.
!! These all resolve to MPP_FILL_INT
  !!Should catch real "default_fill"(MPP_FILL_DOUBLE)
  if (mask_val == MPP_FILL_DOUBLE ) then !this is FMS variable field default fill
     ! we've packed an MPP_FILL_
     imask_val = MPP_FILL_INT
 !Current NETCDF fill values (AKA MPP_FILL_*) designed towards CEILING(MPP_FILL_{FLOAT,DOUBLE},kind=4byte)=MPP_FILL_INT
  else if ( CEILING(mask_val, i4_kind) == MPP_FILL_INT ) then
     ! we've also packed an MPP_FILL_
     imask_val = MPP_FILL_INT
! Secondary Logic:
!! We've done something dangerous
  else
     i8tmp = TRANSFER(mask_val , i8tmp )
     i4tmp = TRANSFER(mask_val , i4tmp )
     r4tmp = TRANSFER(mask_val , r4tmp )
     if ( i8tmp == MPP_FILL_INT ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(i4tmp == MPP_FILL_INT) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(r4tmp == MPP_FILL_DOUBLE) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else
        ! we have no idea what this is
        ! construct detailed errStr
        errStr = "mpp_chksum: mpp_chksum_i"
        write(unit=tmpStr1,fmt="(I1)") KIND(var)
        write(unit=tmpstr2,fmt="(I1)") SIZE(SHAPE(var))
        errStr = errStr // tmpStr1 // "_" // tmpstr2 // "d_rmask passed int var with REAL("
        write(unit=tmpstr3,fmt="(I1)") KIND(mask_val)
        errStr = errStr // tmpstr3 // ") mask_val="
        write(unit=tmpstr4,fmt=*) mask_val
        errStr = errStr // trim(tmpstr4) // "has been called with these strange values. Check your KINDS, "// &
              "_FillValue, pack and mask_val. "// &
              "Hint: Try being explicit and using MPP_FILL_{INT,FLOAT,DOUBLE}. "// &
              "Continuing by using the default MPP_FILL_INT. " // &
              "THIS WILL BE FATAL IN THE FUTURE!"
        call mpp_error(WARNING, trim(errStr) )

        imask_val = MPP_FILL_INT
     end if
  end if

  mpp_chksum_i4_3d_rmask = mpp_chksum(var,pelist,mask_val=imask_val)

  return

end function mpp_chksum_i4_3d_rmask
# 119 "mpp/include/mpp_comm.inc" 2

# 129 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum_int.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> Calculates integer checksum over pelist
function mpp_chksum_i4_4d( var, pelist, mask_val )
  integer(i8_kind) :: mpp_chksum_i4_4d
  integer(i4_kind), intent(in) :: var (:,:,:,:)
  integer, optional :: pelist(:)
  integer(i4_kind), intent(in), optional :: mask_val

  if ( PRESENT(mask_val) ) then
     !PACK on var/=mask_val ignores values in var
     !equiv to setting those values=0, but on sparse arrays
     !pack should return much smaller array to sum
     mpp_chksum_i4_4d = sum( INT( PACK(var,var/=mask_val),i8_kind) )
  else
     mpp_chksum_i4_4d = sum(INT(var,i8_kind))
  end if

      call mpp_sum( mpp_chksum_i4_4d, pelist )
      return

    end function mpp_chksum_i4_4d


!> Handles real mask for easier implementation
function mpp_chksum_i4_4d_rmask( var, pelist, mask_val )
  integer(KIND=i8_kind) :: mpp_chksum_i4_4d_rmask
  integer(i4_kind), intent(in) :: var (:,:,:,:)
  integer, optional :: pelist(:)
  real, intent(in) :: mask_val
  integer(KIND(var))::imask_val
  integer(KIND=i4_kind)::i4tmp(2)=0
  real(KIND=r4_kind)::r4tmp(2)=0
  integer(KIND=i8_kind) :: i8tmp=0
  !high fidelity error message
  character(LEN=1) :: tmpStr1,tmpStr2,tmpStr3
  character(LEN=32) :: tmpStr4
  character(LEN=512) :: errStr

! Primary Logic: These first two are the "expected" branches.
!! These all resolve to MPP_FILL_INT
  !!Should catch real "default_fill"(MPP_FILL_DOUBLE)
  if (mask_val == MPP_FILL_DOUBLE ) then !this is FMS variable field default fill
     ! we've packed an MPP_FILL_
     imask_val = MPP_FILL_INT
 !Current NETCDF fill values (AKA MPP_FILL_*) designed towards CEILING(MPP_FILL_{FLOAT,DOUBLE},kind=4byte)=MPP_FILL_INT
  else if ( CEILING(mask_val, i4_kind) == MPP_FILL_INT ) then
     ! we've also packed an MPP_FILL_
     imask_val = MPP_FILL_INT
! Secondary Logic:
!! We've done something dangerous
  else
     i8tmp = TRANSFER(mask_val , i8tmp )
     i4tmp = TRANSFER(mask_val , i4tmp )
     r4tmp = TRANSFER(mask_val , r4tmp )
     if ( i8tmp == MPP_FILL_INT ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(i4tmp == MPP_FILL_INT) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(r4tmp == MPP_FILL_DOUBLE) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else
        ! we have no idea what this is
        ! construct detailed errStr
        errStr = "mpp_chksum: mpp_chksum_i"
        write(unit=tmpStr1,fmt="(I1)") KIND(var)
        write(unit=tmpstr2,fmt="(I1)") SIZE(SHAPE(var))
        errStr = errStr // tmpStr1 // "_" // tmpstr2 // "d_rmask passed int var with REAL("
        write(unit=tmpstr3,fmt="(I1)") KIND(mask_val)
        errStr = errStr // tmpstr3 // ") mask_val="
        write(unit=tmpstr4,fmt=*) mask_val
        errStr = errStr // trim(tmpstr4) // "has been called with these strange values. Check your KINDS, "// &
              "_FillValue, pack and mask_val. "// &
              "Hint: Try being explicit and using MPP_FILL_{INT,FLOAT,DOUBLE}. "// &
              "Continuing by using the default MPP_FILL_INT. " // &
              "THIS WILL BE FATAL IN THE FUTURE!"
        call mpp_error(WARNING, trim(errStr) )

        imask_val = MPP_FILL_INT
     end if
  end if

  mpp_chksum_i4_4d_rmask = mpp_chksum(var,pelist,mask_val=imask_val)

  return

end function mpp_chksum_i4_4d_rmask
# 129 "mpp/include/mpp_comm.inc" 2

# 139 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum_int.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> Calculates integer checksum over pelist
function mpp_chksum_i4_5d( var, pelist, mask_val )
  integer(i8_kind) :: mpp_chksum_i4_5d
  integer(i4_kind), intent(in) :: var (:,:,:,:,:)
  integer, optional :: pelist(:)
  integer(i4_kind), intent(in), optional :: mask_val

  if ( PRESENT(mask_val) ) then
     !PACK on var/=mask_val ignores values in var
     !equiv to setting those values=0, but on sparse arrays
     !pack should return much smaller array to sum
     mpp_chksum_i4_5d = sum( INT( PACK(var,var/=mask_val),i8_kind) )
  else
     mpp_chksum_i4_5d = sum(INT(var,i8_kind))
  end if

      call mpp_sum( mpp_chksum_i4_5d, pelist )
      return

    end function mpp_chksum_i4_5d


!> Handles real mask for easier implementation
function mpp_chksum_i4_5d_rmask( var, pelist, mask_val )
  integer(KIND=i8_kind) :: mpp_chksum_i4_5d_rmask
  integer(i4_kind), intent(in) :: var (:,:,:,:,:)
  integer, optional :: pelist(:)
  real, intent(in) :: mask_val
  integer(KIND(var))::imask_val
  integer(KIND=i4_kind)::i4tmp(2)=0
  real(KIND=r4_kind)::r4tmp(2)=0
  integer(KIND=i8_kind) :: i8tmp=0
  !high fidelity error message
  character(LEN=1) :: tmpStr1,tmpStr2,tmpStr3
  character(LEN=32) :: tmpStr4
  character(LEN=512) :: errStr

! Primary Logic: These first two are the "expected" branches.
!! These all resolve to MPP_FILL_INT
  !!Should catch real "default_fill"(MPP_FILL_DOUBLE)
  if (mask_val == MPP_FILL_DOUBLE ) then !this is FMS variable field default fill
     ! we've packed an MPP_FILL_
     imask_val = MPP_FILL_INT
 !Current NETCDF fill values (AKA MPP_FILL_*) designed towards CEILING(MPP_FILL_{FLOAT,DOUBLE},kind=4byte)=MPP_FILL_INT
  else if ( CEILING(mask_val, i4_kind) == MPP_FILL_INT ) then
     ! we've also packed an MPP_FILL_
     imask_val = MPP_FILL_INT
! Secondary Logic:
!! We've done something dangerous
  else
     i8tmp = TRANSFER(mask_val , i8tmp )
     i4tmp = TRANSFER(mask_val , i4tmp )
     r4tmp = TRANSFER(mask_val , r4tmp )
     if ( i8tmp == MPP_FILL_INT ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(i4tmp == MPP_FILL_INT) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else if ( ANY(r4tmp == MPP_FILL_DOUBLE) ) then
        ! we've packed an MPP_FILL_
        imask_val = MPP_FILL_INT
     else
        ! we have no idea what this is
        ! construct detailed errStr
        errStr = "mpp_chksum: mpp_chksum_i"
        write(unit=tmpStr1,fmt="(I1)") KIND(var)
        write(unit=tmpstr2,fmt="(I1)") SIZE(SHAPE(var))
        errStr = errStr // tmpStr1 // "_" // tmpstr2 // "d_rmask passed int var with REAL("
        write(unit=tmpstr3,fmt="(I1)") KIND(mask_val)
        errStr = errStr // tmpstr3 // ") mask_val="
        write(unit=tmpstr4,fmt=*) mask_val
        errStr = errStr // trim(tmpstr4) // "has been called with these strange values. Check your KINDS, "// &
              "_FillValue, pack and mask_val. "// &
              "Hint: Try being explicit and using MPP_FILL_{INT,FLOAT,DOUBLE}. "// &
              "Continuing by using the default MPP_FILL_INT. " // &
              "THIS WILL BE FATAL IN THE FUTURE!"
        call mpp_error(WARNING, trim(errStr) )

        imask_val = MPP_FILL_INT
     end if
  end if

  mpp_chksum_i4_5d_rmask = mpp_chksum(var,pelist,mask_val=imask_val)

  return

end function mpp_chksum_i4_5d_rmask
# 139 "mpp/include/mpp_comm.inc" 2








# 1 "mpp/include/mpp_chksum_scalar.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> @brief Wrapper routine for scalar checksums
!!
!> mold is a dummy array to be used by TRANSFER()
!! must be same TYPE as result
!! result is i8_kind, which will actually be int ifdef no_8byte_integers
!! mold and mask_val must be same numBytes, otherwise undefined behavior
function mpp_chksum_r8_0d( var, pelist, mask_val )
      integer(i8_kind) :: mpp_chksum_r8_0d
      real(r8_kind), intent(in) :: var
      integer, intent(in), optional :: pelist(:)
      integer(i8_kind) :: mold(1)
  real(r8_kind), intent(in), optional :: mask_val
      pointer( p, mold )

      p = LOC(var)

  if ( PRESENT(mask_val) ) then
     mpp_chksum_r8_0d = mpp_chksum( mold, pelist, TRANSFER(mask_val, mold(1)) )
  else
      mpp_chksum_r8_0d = mpp_chksum( mold, pelist )
  end if
      return
    end function mpp_chksum_r8_0d
# 147 "mpp/include/mpp_comm.inc" 2

# 157 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> Wrapper routine for @ref mpp_chksum interface
!!
!> @returns i8_kind checksum of var, which will actually be int if no_8byte_integers is defined
function mpp_chksum_r8_1d( var, pelist , mask_val)
  integer(i8_kind) :: mpp_chksum_r8_1d
  integer(i8_kind) :: mold(1) !< Mold is a dummy array to be used by TRANSFER(),
                                         !! must be same TYPE as result
  real(r8_kind), intent(in) :: var (:)
  integer, intent(in), optional :: pelist(:)
  real(r8_kind), intent(in),optional :: mask_val !< optional mask_val is masked away in checksum_int.h
                                             !! function via PACK()

  if ( PRESENT(mask_val) ) then
     mpp_chksum_r8_1d = mpp_chksum( TRANSFER(var,mold), pelist, &
          mask_val= TRANSFER(mask_val,mold(1) ) )
  else
      mpp_chksum_r8_1d = mpp_chksum( TRANSFER(var,mold), pelist )
  end if

  return
end function mpp_chksum_r8_1d
# 157 "mpp/include/mpp_comm.inc" 2

# 167 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> Wrapper routine for @ref mpp_chksum interface
!!
!> @returns i8_kind checksum of var, which will actually be int if no_8byte_integers is defined
function mpp_chksum_r8_2d( var, pelist , mask_val)
  integer(i8_kind) :: mpp_chksum_r8_2d
  integer(i8_kind) :: mold(1) !< Mold is a dummy array to be used by TRANSFER(),
                                         !! must be same TYPE as result
  real(r8_kind), intent(in) :: var (:,:)
  integer, intent(in), optional :: pelist(:)
  real(r8_kind), intent(in),optional :: mask_val !< optional mask_val is masked away in checksum_int.h
                                             !! function via PACK()

  if ( PRESENT(mask_val) ) then
     mpp_chksum_r8_2d = mpp_chksum( TRANSFER(var,mold), pelist, &
          mask_val= TRANSFER(mask_val,mold(1) ) )
  else
      mpp_chksum_r8_2d = mpp_chksum( TRANSFER(var,mold), pelist )
  end if

  return
end function mpp_chksum_r8_2d
# 167 "mpp/include/mpp_comm.inc" 2

# 177 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> Wrapper routine for @ref mpp_chksum interface
!!
!> @returns i8_kind checksum of var, which will actually be int if no_8byte_integers is defined
function mpp_chksum_r8_3d( var, pelist , mask_val)
  integer(i8_kind) :: mpp_chksum_r8_3d
  integer(i8_kind) :: mold(1) !< Mold is a dummy array to be used by TRANSFER(),
                                         !! must be same TYPE as result
  real(r8_kind), intent(in) :: var (:,:,:)
  integer, intent(in), optional :: pelist(:)
  real(r8_kind), intent(in),optional :: mask_val !< optional mask_val is masked away in checksum_int.h
                                             !! function via PACK()

  if ( PRESENT(mask_val) ) then
     mpp_chksum_r8_3d = mpp_chksum( TRANSFER(var,mold), pelist, &
          mask_val= TRANSFER(mask_val,mold(1) ) )
  else
      mpp_chksum_r8_3d = mpp_chksum( TRANSFER(var,mold), pelist )
  end if

  return
end function mpp_chksum_r8_3d
# 177 "mpp/include/mpp_comm.inc" 2

# 187 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> Wrapper routine for @ref mpp_chksum interface
!!
!> @returns i8_kind checksum of var, which will actually be int if no_8byte_integers is defined
function mpp_chksum_r8_4d( var, pelist , mask_val)
  integer(i8_kind) :: mpp_chksum_r8_4d
  integer(i8_kind) :: mold(1) !< Mold is a dummy array to be used by TRANSFER(),
                                         !! must be same TYPE as result
  real(r8_kind), intent(in) :: var (:,:,:,:)
  integer, intent(in), optional :: pelist(:)
  real(r8_kind), intent(in),optional :: mask_val !< optional mask_val is masked away in checksum_int.h
                                             !! function via PACK()

  if ( PRESENT(mask_val) ) then
     mpp_chksum_r8_4d = mpp_chksum( TRANSFER(var,mold), pelist, &
          mask_val= TRANSFER(mask_val,mold(1) ) )
  else
      mpp_chksum_r8_4d = mpp_chksum( TRANSFER(var,mold), pelist )
  end if

  return
end function mpp_chksum_r8_4d
# 187 "mpp/include/mpp_comm.inc" 2

# 197 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> Wrapper routine for @ref mpp_chksum interface
!!
!> @returns i8_kind checksum of var, which will actually be int if no_8byte_integers is defined
function mpp_chksum_r8_5d( var, pelist , mask_val)
  integer(i8_kind) :: mpp_chksum_r8_5d
  integer(i8_kind) :: mold(1) !< Mold is a dummy array to be used by TRANSFER(),
                                         !! must be same TYPE as result
  real(r8_kind), intent(in) :: var (:,:,:,:,:)
  integer, intent(in), optional :: pelist(:)
  real(r8_kind), intent(in),optional :: mask_val !< optional mask_val is masked away in checksum_int.h
                                             !! function via PACK()

  if ( PRESENT(mask_val) ) then
     mpp_chksum_r8_5d = mpp_chksum( TRANSFER(var,mold), pelist, &
          mask_val= TRANSFER(mask_val,mold(1) ) )
  else
      mpp_chksum_r8_5d = mpp_chksum( TRANSFER(var,mold), pelist )
  end if

  return
end function mpp_chksum_r8_5d
# 197 "mpp/include/mpp_comm.inc" 2

# 259 "mpp/include/mpp_comm.inc"

# 269 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum_scalar.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> @brief Wrapper routine for scalar checksums
!!
!> mold is a dummy array to be used by TRANSFER()
!! must be same TYPE as result
!! result is i8_kind, which will actually be int ifdef no_8byte_integers
!! mold and mask_val must be same numBytes, otherwise undefined behavior
function mpp_chksum_r4_0d( var, pelist, mask_val )
      integer(i8_kind) :: mpp_chksum_r4_0d
      real(r4_kind), intent(in) :: var
      integer, intent(in), optional :: pelist(:)
      integer(i8_kind) :: mold(1)
  real(r4_kind), intent(in), optional :: mask_val
      pointer( p, mold )

      p = LOC(var)

  if ( PRESENT(mask_val) ) then
     mpp_chksum_r4_0d = mpp_chksum( mold, pelist, TRANSFER(mask_val, mold(1)) )
  else
      mpp_chksum_r4_0d = mpp_chksum( mold, pelist )
  end if
      return
    end function mpp_chksum_r4_0d
# 269 "mpp/include/mpp_comm.inc" 2

# 279 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> Wrapper routine for @ref mpp_chksum interface
!!
!> @returns i8_kind checksum of var, which will actually be int if no_8byte_integers is defined
function mpp_chksum_r4_1d( var, pelist , mask_val)
  integer(i8_kind) :: mpp_chksum_r4_1d
  integer(i4_kind) :: mold(1) !< Mold is a dummy array to be used by TRANSFER(),
                                         !! must be same TYPE as result
  real(r4_kind), intent(in) :: var (:)
  integer, intent(in), optional :: pelist(:)
  real(r4_kind), intent(in),optional :: mask_val !< optional mask_val is masked away in checksum_int.h
                                             !! function via PACK()

  if ( PRESENT(mask_val) ) then
     mpp_chksum_r4_1d = mpp_chksum( TRANSFER(var,mold), pelist, &
          mask_val= TRANSFER(mask_val,mold(1) ) )
  else
      mpp_chksum_r4_1d = mpp_chksum( TRANSFER(var,mold), pelist )
  end if

  return
end function mpp_chksum_r4_1d
# 279 "mpp/include/mpp_comm.inc" 2

# 289 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> Wrapper routine for @ref mpp_chksum interface
!!
!> @returns i8_kind checksum of var, which will actually be int if no_8byte_integers is defined
function mpp_chksum_r4_2d( var, pelist , mask_val)
  integer(i8_kind) :: mpp_chksum_r4_2d
  integer(i4_kind) :: mold(1) !< Mold is a dummy array to be used by TRANSFER(),
                                         !! must be same TYPE as result
  real(r4_kind), intent(in) :: var (:,:)
  integer, intent(in), optional :: pelist(:)
  real(r4_kind), intent(in),optional :: mask_val !< optional mask_val is masked away in checksum_int.h
                                             !! function via PACK()

  if ( PRESENT(mask_val) ) then
     mpp_chksum_r4_2d = mpp_chksum( TRANSFER(var,mold), pelist, &
          mask_val= TRANSFER(mask_val,mold(1) ) )
  else
      mpp_chksum_r4_2d = mpp_chksum( TRANSFER(var,mold), pelist )
  end if

  return
end function mpp_chksum_r4_2d
# 289 "mpp/include/mpp_comm.inc" 2

# 299 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> Wrapper routine for @ref mpp_chksum interface
!!
!> @returns i8_kind checksum of var, which will actually be int if no_8byte_integers is defined
function mpp_chksum_r4_3d( var, pelist , mask_val)
  integer(i8_kind) :: mpp_chksum_r4_3d
  integer(i4_kind) :: mold(1) !< Mold is a dummy array to be used by TRANSFER(),
                                         !! must be same TYPE as result
  real(r4_kind), intent(in) :: var (:,:,:)
  integer, intent(in), optional :: pelist(:)
  real(r4_kind), intent(in),optional :: mask_val !< optional mask_val is masked away in checksum_int.h
                                             !! function via PACK()

  if ( PRESENT(mask_val) ) then
     mpp_chksum_r4_3d = mpp_chksum( TRANSFER(var,mold), pelist, &
          mask_val= TRANSFER(mask_val,mold(1) ) )
  else
      mpp_chksum_r4_3d = mpp_chksum( TRANSFER(var,mold), pelist )
  end if

  return
end function mpp_chksum_r4_3d
# 299 "mpp/include/mpp_comm.inc" 2

# 309 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> Wrapper routine for @ref mpp_chksum interface
!!
!> @returns i8_kind checksum of var, which will actually be int if no_8byte_integers is defined
function mpp_chksum_r4_4d( var, pelist , mask_val)
  integer(i8_kind) :: mpp_chksum_r4_4d
  integer(i4_kind) :: mold(1) !< Mold is a dummy array to be used by TRANSFER(),
                                         !! must be same TYPE as result
  real(r4_kind), intent(in) :: var (:,:,:,:)
  integer, intent(in), optional :: pelist(:)
  real(r4_kind), intent(in),optional :: mask_val !< optional mask_val is masked away in checksum_int.h
                                             !! function via PACK()

  if ( PRESENT(mask_val) ) then
     mpp_chksum_r4_4d = mpp_chksum( TRANSFER(var,mold), pelist, &
          mask_val= TRANSFER(mask_val,mold(1) ) )
  else
      mpp_chksum_r4_4d = mpp_chksum( TRANSFER(var,mold), pelist )
  end if

  return
end function mpp_chksum_r4_4d
# 309 "mpp/include/mpp_comm.inc" 2

# 319 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_chksum.fh" 1
! -*-f90-*-

!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> Wrapper routine for @ref mpp_chksum interface
!!
!> @returns i8_kind checksum of var, which will actually be int if no_8byte_integers is defined
function mpp_chksum_r4_5d( var, pelist , mask_val)
  integer(i8_kind) :: mpp_chksum_r4_5d
  integer(i4_kind) :: mold(1) !< Mold is a dummy array to be used by TRANSFER(),
                                         !! must be same TYPE as result
  real(r4_kind), intent(in) :: var (:,:,:,:,:)
  integer, intent(in), optional :: pelist(:)
  real(r4_kind), intent(in),optional :: mask_val !< optional mask_val is masked away in checksum_int.h
                                             !! function via PACK()

  if ( PRESENT(mask_val) ) then
     mpp_chksum_r4_5d = mpp_chksum( TRANSFER(var,mold), pelist, &
          mask_val= TRANSFER(mask_val,mold(1) ) )
  else
      mpp_chksum_r4_5d = mpp_chksum( TRANSFER(var,mold), pelist )
  end if

  return
end function mpp_chksum_r4_5d
# 319 "mpp/include/mpp_comm.inc" 2

# 381 "mpp/include/mpp_comm.inc"

!#################################################
# 400 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_gather.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_gather_logical_1d(sbuf, rbuf, pelist)
! JWD: Did not create mpp_gather_2d because have no requirement for it
! JWD: See mpp_gather_2dv below
   logical, dimension(:),    intent(in) :: sbuf
   logical, dimension(:), intent(inout) :: rbuf
   integer,   dimension(:),    intent(in), optional :: pelist(:)

   integer :: cnt, l, nproc, op_root, ierr
   integer, allocatable :: pelist2(:)

   if( .NOT.module_is_initialized ) call mpp_error( FATAL, 'MPP_GATHER_1D_: You must first call mpp_init.' )

!  If pelist is provided, the first position must be the operation root, w.r.t. new comm, op_root = 0
   if(PRESENT(pelist))then
      if(.not.ANY(mpp_pe().eq.pelist(:))) return
      nproc = size(pelist)
      allocate(pelist2(nproc))
      pelist2 = pelist
   else
      nproc = mpp_npes()
      allocate(pelist2(nproc))
      pelist2 = (/ (l, l=root_pe, nproc-1+root_pe) /)
   endif
   op_root = pelist2(1)

   cnt = size(sbuf(:))
   if((mpp_pe().eq.op_root).AND.(size(rbuf(:)) < cnt*nproc)) call mpp_error(FATAL, &
          "MPP_GATHER_1D_: size(rbuf) must be at least npes*size(sbuf) ")

   call mpp_gather( sbuf, rbuf, size(sbuf), op_root, pelist2, ierr )

   call mpp_sync_self()
   deallocate(pelist2)
end subroutine mpp_gather_logical_1d

subroutine mpp_gather_logical_1dv(sbuf, ssize, rbuf, rsize, pelist)
   logical, dimension(:),    intent(in) :: sbuf
   logical, dimension(:), intent(inout) :: rbuf
   integer,                    intent(in) :: ssize
   integer,   dimension(:),    intent(in) :: rsize
   integer,   dimension(:),    intent(in), optional :: pelist(:)

   integer :: l, nproc, op_root, ierr
   integer, dimension(:), allocatable :: displs
   integer, dimension(:), allocatable :: pelist2

!  If pelist is provided, the first position must be
!  the operation root
   if(PRESENT(pelist))then
      nproc = size(pelist)
      allocate(pelist2(nproc))
      pelist2 = pelist
   else
      nproc = mpp_npes()
      allocate(pelist2(nproc))
      pelist2 = (/ (l, l=root_pe, nproc-1+root_pe) /)
   endif
   op_root = pelist2(1)

   if(pe .eq. op_root) then
      allocate(displs(nproc))

      displs(1) = 0
      do l = 2, nproc
         displs(l) = displs(l-1) + rsize(l-1)
      enddo
   else
      allocate(displs(1))
   endif

   call mpp_gather( sbuf, ssize, rbuf, rsize, displs, op_root, pelist2, ierr )

   call mpp_sync_self()
   deallocate(pelist2)
   deallocate(displs)
end subroutine mpp_gather_logical_1dv

subroutine mpp_gather_pelist_logical_2d(is, ie, js, je, pelist, array_seg, gather_data, is_root_pe, &
                                 ishift, jshift)
   integer,                                       intent(in)    :: is, ie, js, je
   integer,   dimension(:),                       intent(in)    :: pelist
   logical, dimension(is:ie,js:je),     target, intent(in)    :: array_seg
   logical, dimension(:,:), contiguous, target, intent(inout) :: gather_data
   logical,                                       intent(in)    :: is_root_pe
   integer,   optional,                           intent(in)    :: ishift, jshift

   integer, dimension(2)  :: axis_to_storage

   axis_to_storage = (/1,2/)

   call mpp_gather(is, ie, js, je, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                   ishift, jshift)
   return

end subroutine mpp_gather_pelist_logical_2d

subroutine mpp_gather_pelist_logical_gen_2d(is, ie, js, je, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                                     ishift, jshift)
   integer,                                       intent(in)    :: is, ie, js, je
   integer,   dimension(:),                       intent(in)    :: pelist
   logical, dimension(:,:), contiguous, target, intent(in)    :: array_seg
   logical, dimension(:,:), contiguous, target, intent(inout) :: gather_data
   integer,   dimension(2),                       intent(in)    :: axis_to_storage
   logical,                                       intent(in)    :: is_root_pe
   integer,   optional,                           intent(in)    :: ishift, jshift

   logical, pointer ::  arr3D(:,:,:)
   logical, pointer :: data3D(:,:,:)

   arr3D(1:size(array_seg,1),1:size(array_seg,2),1:1) => array_seg
   if (is_root_pe) then
     data3D(1:size(gather_data,1),1:size(gather_data,2),1:1) => gather_data
   else
     data3D => null()
   endif

   call mpp_gather(is, ie, js, je, 1, pelist, arr3D, data3D, [axis_to_storage, 3], is_root_pe, &
                   ishift, jshift)
   return

end subroutine mpp_gather_pelist_logical_gen_2d


subroutine mpp_gather_pelist_logical_3d(is, ie, js, je, nk, pelist, array_seg, gather_data, is_root_pe, &
                                 ishift, jshift)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   logical, dimension(is:ie,js:je,1:nk), intent(in)    :: array_seg
   logical, dimension(:,:,:),            intent(inout) :: gather_data
   logical,                                intent(in)    :: is_root_pe
   integer,   optional,                    intent(in)    :: ishift, jshift

   integer, dimension(3) :: axis_to_storage

   axis_to_storage = (/1, 2, 3/)

   call mpp_gather(is, ie, js, je, nk, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                   ishift, jshift)
   return

end subroutine mpp_gather_pelist_logical_3d

subroutine mpp_gather_pelist_logical_gen_3d(is, ie, js, je, nk, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                                     ishift, jshift)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   logical, dimension(:,:,:),            intent(in)    :: array_seg
   logical, dimension(:,:,:),            intent(inout) :: gather_data
   integer,   dimension(3),                intent(in)    :: axis_to_storage
   logical,                                intent(in)    :: is_root_pe
   integer,   optional,                    intent(in)    :: ishift, jshift

   integer :: root_pe, root_pe_test
   integer :: k, us, ue, vs, ve, ws, we
   integer :: i1, i2, j1, j2, ioff, joff
   integer :: base_idx, send_count, msg_start
   integer :: blocksize_u, blocksize_v, blocksize_w, blocksize
   integer, dimension(3) :: start_idx, stop_idx, storage_to_axis
   integer, dimension(:), allocatable :: gind, counts
   logical, dimension(:), allocatable :: rbuf

   if (.not.ANY(mpp_pe().eq.pelist(:))) return

   ! Check axis_to_storage is a permutation of 1..3
   if ( any(axis_to_storage < 1) .or. any(axis_to_storage > 3) ) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): axis_to_storage entries must be in {1,2,3}")
   if ( axis_to_storage(1) == axis_to_storage(2) .or. axis_to_storage(1) == axis_to_storage(3) &
        .or. axis_to_storage(2) == axis_to_storage(3) ) &
     call mpp_error(FATAL, "fms_io(mpp_gather_pelist): axis_to_storage must be a permutation of 1,2,3")

   if (is_root_pe) then
     root_pe = mpp_pe()
     root_pe_test = 999
     if (.not.ANY(pelist(:).eq.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): root_pe not a member of pelist")
   else
     root_pe = 0
     root_pe_test = -999
   endif
! need this check in case MPI-rank 0 is a member of the pelist
   call mpp_max(root_pe_test, pelist)
   if (root_pe_test.lt.0) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): root_pe not specified or not a member of the pelist")
! need to make sure only one root_pe has been specified
   call mpp_sum(root_pe, pelist)
   if ((is_root_pe) .and. (mpp_pe().ne.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): too many root_pes specified")

   ioff=0
   joff=0
   if (present(ishift)) ioff=ishift
   if (present(jshift)) joff=jshift

   ! Initialize storage dim to axis map
   storage_to_axis(axis_to_storage(1)) = 1
   storage_to_axis(axis_to_storage(2)) = 2
   storage_to_axis(axis_to_storage(3)) = 3

   ! gather indices into global index on root_pe
   if (is_root_pe) then
     allocate(gind(4*size(pelist)))
   else
     allocate(gind(1))
   endif
   call mpp_gather((/is, ie, js, je/), gind, pelist)

   ! Compute recv counts and allocate 1d recv buffer (rbuf)
   if (is_root_pe) then
      allocate(counts(size(pelist)))

      do k = 1, size(pelist)
         base_idx = 4*(k-1)
         i1 = gind( base_idx + 1 ) ;; i2 = gind( base_idx + 2 )
         j1 = gind( base_idx + 3 ) ;; j2 = gind( base_idx + 4 )
         counts(k) = (i2 - i1 + 1) * (j2 - j1 + 1) * nk
      enddo

      allocate(rbuf(sum(counts)))
   else
      ! Non-root: MPI ignores recv args, but they must still be valid actual arguements
      allocate(counts(1))
      counts = 0
      allocate(rbuf(1))
   endif

   send_count = (ie-is+1)*(je-js+1)*nk

   ! Get generalized stop indicies for array_seg
   stop_idx = (/ie-is+1, je-js+1, nk/)
   ue = stop_idx(storage_to_axis(1))
   ve = stop_idx(storage_to_axis(2))
   we = stop_idx(storage_to_axis(3))

   ! gather data into 1d recv buffer
   call mpp_gather(reshape(array_seg(1:ue,1:ve,1:we),[send_count]), send_count, rbuf, counts, pelist)

   ! Unpack recv buffer into return array (gather_data)
   if (is_root_pe) then
      msg_start = 1
      do k = 1, size(pelist)
         base_idx = 4*(k-1)
         i1 = gind( base_idx + 1 ) + ioff ;; i2 = gind( base_idx + 2 ) + ioff
         j1 = gind( base_idx + 3 ) + joff ;; j2 = gind( base_idx + 4 ) + joff

         ! Get generalized start/stop indicies
         start_idx = (/i1,j1,1/)
         stop_idx  = (/i2,j2,nk/)

         us = start_idx(storage_to_axis(1)) ;; ue = stop_idx(storage_to_axis(1))
         vs = start_idx(storage_to_axis(2)) ;; ve = stop_idx(storage_to_axis(2))
         ws = start_idx(storage_to_axis(3)) ;; we = stop_idx(storage_to_axis(3))

         ! Compute block sizes
         blocksize_u = ue - us + 1
         blocksize_v = ve - vs + 1
         blocksize_w = we - ws + 1
         blocksize   = blocksize_u * blocksize_v * blocksize_w

         gather_data(us:ue, vs:ve, ws:we) = reshape(rbuf(msg_start:msg_start+blocksize-1), &
                                                         [blocksize_u, blocksize_v, blocksize_w])

         msg_start = msg_start + blocksize
      enddo

      deallocate(gind)
   endif

   deallocate(rbuf, counts)

   call mpp_sync_self()

end subroutine mpp_gather_pelist_logical_gen_3d
# 400 "mpp/include/mpp_comm.inc" 2

# 418 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_gather.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_gather_int4_1d(sbuf, rbuf, pelist)
! JWD: Did not create mpp_gather_2d because have no requirement for it
! JWD: See mpp_gather_2dv below
   integer(i4_kind), dimension(:),    intent(in) :: sbuf
   integer(i4_kind), dimension(:), intent(inout) :: rbuf
   integer,   dimension(:),    intent(in), optional :: pelist(:)

   integer :: cnt, l, nproc, op_root, ierr
   integer, allocatable :: pelist2(:)

   if( .NOT.module_is_initialized ) call mpp_error( FATAL, 'MPP_GATHER_1D_: You must first call mpp_init.' )

!  If pelist is provided, the first position must be the operation root, w.r.t. new comm, op_root = 0
   if(PRESENT(pelist))then
      if(.not.ANY(mpp_pe().eq.pelist(:))) return
      nproc = size(pelist)
      allocate(pelist2(nproc))
      pelist2 = pelist
   else
      nproc = mpp_npes()
      allocate(pelist2(nproc))
      pelist2 = (/ (l, l=root_pe, nproc-1+root_pe) /)
   endif
   op_root = pelist2(1)

   cnt = size(sbuf(:))
   if((mpp_pe().eq.op_root).AND.(size(rbuf(:)) < cnt*nproc)) call mpp_error(FATAL, &
          "MPP_GATHER_1D_: size(rbuf) must be at least npes*size(sbuf) ")

   call mpp_gather( sbuf, rbuf, size(sbuf), op_root, pelist2, ierr )

   call mpp_sync_self()
   deallocate(pelist2)
end subroutine mpp_gather_int4_1d

subroutine mpp_gather_int4_1dv(sbuf, ssize, rbuf, rsize, pelist)
   integer(i4_kind), dimension(:),    intent(in) :: sbuf
   integer(i4_kind), dimension(:), intent(inout) :: rbuf
   integer,                    intent(in) :: ssize
   integer,   dimension(:),    intent(in) :: rsize
   integer,   dimension(:),    intent(in), optional :: pelist(:)

   integer :: l, nproc, op_root, ierr
   integer, dimension(:), allocatable :: displs
   integer, dimension(:), allocatable :: pelist2

!  If pelist is provided, the first position must be
!  the operation root
   if(PRESENT(pelist))then
      nproc = size(pelist)
      allocate(pelist2(nproc))
      pelist2 = pelist
   else
      nproc = mpp_npes()
      allocate(pelist2(nproc))
      pelist2 = (/ (l, l=root_pe, nproc-1+root_pe) /)
   endif
   op_root = pelist2(1)

   if(pe .eq. op_root) then
      allocate(displs(nproc))

      displs(1) = 0
      do l = 2, nproc
         displs(l) = displs(l-1) + rsize(l-1)
      enddo
   else
      allocate(displs(1))
   endif

   call mpp_gather( sbuf, ssize, rbuf, rsize, displs, op_root, pelist2, ierr )

   call mpp_sync_self()
   deallocate(pelist2)
   deallocate(displs)
end subroutine mpp_gather_int4_1dv

subroutine mpp_gather_pelist_int4_2d(is, ie, js, je, pelist, array_seg, gather_data, is_root_pe, &
                                 ishift, jshift)
   integer,                                       intent(in)    :: is, ie, js, je
   integer,   dimension(:),                       intent(in)    :: pelist
   integer(i4_kind), dimension(is:ie,js:je),     target, intent(in)    :: array_seg
   integer(i4_kind), dimension(:,:), contiguous, target, intent(inout) :: gather_data
   logical,                                       intent(in)    :: is_root_pe
   integer,   optional,                           intent(in)    :: ishift, jshift

   integer, dimension(2)  :: axis_to_storage

   axis_to_storage = (/1,2/)

   call mpp_gather(is, ie, js, je, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                   ishift, jshift)
   return

end subroutine mpp_gather_pelist_int4_2d

subroutine mpp_gather_pelist_int4_gen_2d(is, ie, js, je, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                                     ishift, jshift)
   integer,                                       intent(in)    :: is, ie, js, je
   integer,   dimension(:),                       intent(in)    :: pelist
   integer(i4_kind), dimension(:,:), contiguous, target, intent(in)    :: array_seg
   integer(i4_kind), dimension(:,:), contiguous, target, intent(inout) :: gather_data
   integer,   dimension(2),                       intent(in)    :: axis_to_storage
   logical,                                       intent(in)    :: is_root_pe
   integer,   optional,                           intent(in)    :: ishift, jshift

   integer(i4_kind), pointer ::  arr3D(:,:,:)
   integer(i4_kind), pointer :: data3D(:,:,:)

   arr3D(1:size(array_seg,1),1:size(array_seg,2),1:1) => array_seg
   if (is_root_pe) then
     data3D(1:size(gather_data,1),1:size(gather_data,2),1:1) => gather_data
   else
     data3D => null()
   endif

   call mpp_gather(is, ie, js, je, 1, pelist, arr3D, data3D, [axis_to_storage, 3], is_root_pe, &
                   ishift, jshift)
   return

end subroutine mpp_gather_pelist_int4_gen_2d


subroutine mpp_gather_pelist_int4_3d(is, ie, js, je, nk, pelist, array_seg, gather_data, is_root_pe, &
                                 ishift, jshift)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   integer(i4_kind), dimension(is:ie,js:je,1:nk), intent(in)    :: array_seg
   integer(i4_kind), dimension(:,:,:),            intent(inout) :: gather_data
   logical,                                intent(in)    :: is_root_pe
   integer,   optional,                    intent(in)    :: ishift, jshift

   integer, dimension(3) :: axis_to_storage

   axis_to_storage = (/1, 2, 3/)

   call mpp_gather(is, ie, js, je, nk, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                   ishift, jshift)
   return

end subroutine mpp_gather_pelist_int4_3d

subroutine mpp_gather_pelist_int4_gen_3d(is, ie, js, je, nk, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                                     ishift, jshift)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   integer(i4_kind), dimension(:,:,:),            intent(in)    :: array_seg
   integer(i4_kind), dimension(:,:,:),            intent(inout) :: gather_data
   integer,   dimension(3),                intent(in)    :: axis_to_storage
   logical,                                intent(in)    :: is_root_pe
   integer,   optional,                    intent(in)    :: ishift, jshift

   integer :: root_pe, root_pe_test
   integer :: k, us, ue, vs, ve, ws, we
   integer :: i1, i2, j1, j2, ioff, joff
   integer :: base_idx, send_count, msg_start
   integer :: blocksize_u, blocksize_v, blocksize_w, blocksize
   integer, dimension(3) :: start_idx, stop_idx, storage_to_axis
   integer, dimension(:), allocatable :: gind, counts
   integer(i4_kind), dimension(:), allocatable :: rbuf

   if (.not.ANY(mpp_pe().eq.pelist(:))) return

   ! Check axis_to_storage is a permutation of 1..3
   if ( any(axis_to_storage < 1) .or. any(axis_to_storage > 3) ) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): axis_to_storage entries must be in {1,2,3}")
   if ( axis_to_storage(1) == axis_to_storage(2) .or. axis_to_storage(1) == axis_to_storage(3) &
        .or. axis_to_storage(2) == axis_to_storage(3) ) &
     call mpp_error(FATAL, "fms_io(mpp_gather_pelist): axis_to_storage must be a permutation of 1,2,3")

   if (is_root_pe) then
     root_pe = mpp_pe()
     root_pe_test = 999
     if (.not.ANY(pelist(:).eq.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): root_pe not a member of pelist")
   else
     root_pe = 0
     root_pe_test = -999
   endif
! need this check in case MPI-rank 0 is a member of the pelist
   call mpp_max(root_pe_test, pelist)
   if (root_pe_test.lt.0) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): root_pe not specified or not a member of the pelist")
! need to make sure only one root_pe has been specified
   call mpp_sum(root_pe, pelist)
   if ((is_root_pe) .and. (mpp_pe().ne.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): too many root_pes specified")

   ioff=0
   joff=0
   if (present(ishift)) ioff=ishift
   if (present(jshift)) joff=jshift

   ! Initialize storage dim to axis map
   storage_to_axis(axis_to_storage(1)) = 1
   storage_to_axis(axis_to_storage(2)) = 2
   storage_to_axis(axis_to_storage(3)) = 3

   ! gather indices into global index on root_pe
   if (is_root_pe) then
     allocate(gind(4*size(pelist)))
   else
     allocate(gind(1))
   endif
   call mpp_gather((/is, ie, js, je/), gind, pelist)

   ! Compute recv counts and allocate 1d recv buffer (rbuf)
   if (is_root_pe) then
      allocate(counts(size(pelist)))

      do k = 1, size(pelist)
         base_idx = 4*(k-1)
         i1 = gind( base_idx + 1 ) ;; i2 = gind( base_idx + 2 )
         j1 = gind( base_idx + 3 ) ;; j2 = gind( base_idx + 4 )
         counts(k) = (i2 - i1 + 1) * (j2 - j1 + 1) * nk
      enddo

      allocate(rbuf(sum(counts)))
   else
      ! Non-root: MPI ignores recv args, but they must still be valid actual arguements
      allocate(counts(1))
      counts = 0
      allocate(rbuf(1))
   endif

   send_count = (ie-is+1)*(je-js+1)*nk

   ! Get generalized stop indicies for array_seg
   stop_idx = (/ie-is+1, je-js+1, nk/)
   ue = stop_idx(storage_to_axis(1))
   ve = stop_idx(storage_to_axis(2))
   we = stop_idx(storage_to_axis(3))

   ! gather data into 1d recv buffer
   call mpp_gather(reshape(array_seg(1:ue,1:ve,1:we),[send_count]), send_count, rbuf, counts, pelist)

   ! Unpack recv buffer into return array (gather_data)
   if (is_root_pe) then
      msg_start = 1
      do k = 1, size(pelist)
         base_idx = 4*(k-1)
         i1 = gind( base_idx + 1 ) + ioff ;; i2 = gind( base_idx + 2 ) + ioff
         j1 = gind( base_idx + 3 ) + joff ;; j2 = gind( base_idx + 4 ) + joff

         ! Get generalized start/stop indicies
         start_idx = (/i1,j1,1/)
         stop_idx  = (/i2,j2,nk/)

         us = start_idx(storage_to_axis(1)) ;; ue = stop_idx(storage_to_axis(1))
         vs = start_idx(storage_to_axis(2)) ;; ve = stop_idx(storage_to_axis(2))
         ws = start_idx(storage_to_axis(3)) ;; we = stop_idx(storage_to_axis(3))

         ! Compute block sizes
         blocksize_u = ue - us + 1
         blocksize_v = ve - vs + 1
         blocksize_w = we - ws + 1
         blocksize   = blocksize_u * blocksize_v * blocksize_w

         gather_data(us:ue, vs:ve, ws:we) = reshape(rbuf(msg_start:msg_start+blocksize-1), &
                                                         [blocksize_u, blocksize_v, blocksize_w])

         msg_start = msg_start + blocksize
      enddo

      deallocate(gind)
   endif

   deallocate(rbuf, counts)

   call mpp_sync_self()

end subroutine mpp_gather_pelist_int4_gen_3d
# 418 "mpp/include/mpp_comm.inc" 2


# 437 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_gather.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_gather_int8_1d(sbuf, rbuf, pelist)
! JWD: Did not create mpp_gather_2d because have no requirement for it
! JWD: See mpp_gather_2dv below
   integer(i8_kind), dimension(:),    intent(in) :: sbuf
   integer(i8_kind), dimension(:), intent(inout) :: rbuf
   integer,   dimension(:),    intent(in), optional :: pelist(:)

   integer :: cnt, l, nproc, op_root, ierr
   integer, allocatable :: pelist2(:)

   if( .NOT.module_is_initialized ) call mpp_error( FATAL, 'MPP_GATHER_1D_: You must first call mpp_init.' )

!  If pelist is provided, the first position must be the operation root, w.r.t. new comm, op_root = 0
   if(PRESENT(pelist))then
      if(.not.ANY(mpp_pe().eq.pelist(:))) return
      nproc = size(pelist)
      allocate(pelist2(nproc))
      pelist2 = pelist
   else
      nproc = mpp_npes()
      allocate(pelist2(nproc))
      pelist2 = (/ (l, l=root_pe, nproc-1+root_pe) /)
   endif
   op_root = pelist2(1)

   cnt = size(sbuf(:))
   if((mpp_pe().eq.op_root).AND.(size(rbuf(:)) < cnt*nproc)) call mpp_error(FATAL, &
          "MPP_GATHER_1D_: size(rbuf) must be at least npes*size(sbuf) ")

   call mpp_gather( sbuf, rbuf, size(sbuf), op_root, pelist2, ierr )

   call mpp_sync_self()
   deallocate(pelist2)
end subroutine mpp_gather_int8_1d

subroutine mpp_gather_int8_1dv(sbuf, ssize, rbuf, rsize, pelist)
   integer(i8_kind), dimension(:),    intent(in) :: sbuf
   integer(i8_kind), dimension(:), intent(inout) :: rbuf
   integer,                    intent(in) :: ssize
   integer,   dimension(:),    intent(in) :: rsize
   integer,   dimension(:),    intent(in), optional :: pelist(:)

   integer :: l, nproc, op_root, ierr
   integer, dimension(:), allocatable :: displs
   integer, dimension(:), allocatable :: pelist2

!  If pelist is provided, the first position must be
!  the operation root
   if(PRESENT(pelist))then
      nproc = size(pelist)
      allocate(pelist2(nproc))
      pelist2 = pelist
   else
      nproc = mpp_npes()
      allocate(pelist2(nproc))
      pelist2 = (/ (l, l=root_pe, nproc-1+root_pe) /)
   endif
   op_root = pelist2(1)

   if(pe .eq. op_root) then
      allocate(displs(nproc))

      displs(1) = 0
      do l = 2, nproc
         displs(l) = displs(l-1) + rsize(l-1)
      enddo
   else
      allocate(displs(1))
   endif

   call mpp_gather( sbuf, ssize, rbuf, rsize, displs, op_root, pelist2, ierr )

   call mpp_sync_self()
   deallocate(pelist2)
   deallocate(displs)
end subroutine mpp_gather_int8_1dv

subroutine mpp_gather_pelist_int8_2d(is, ie, js, je, pelist, array_seg, gather_data, is_root_pe, &
                                 ishift, jshift)
   integer,                                       intent(in)    :: is, ie, js, je
   integer,   dimension(:),                       intent(in)    :: pelist
   integer(i8_kind), dimension(is:ie,js:je),     target, intent(in)    :: array_seg
   integer(i8_kind), dimension(:,:), contiguous, target, intent(inout) :: gather_data
   logical,                                       intent(in)    :: is_root_pe
   integer,   optional,                           intent(in)    :: ishift, jshift

   integer, dimension(2)  :: axis_to_storage

   axis_to_storage = (/1,2/)

   call mpp_gather(is, ie, js, je, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                   ishift, jshift)
   return

end subroutine mpp_gather_pelist_int8_2d

subroutine mpp_gather_pelist_int8_gen_2d(is, ie, js, je, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                                     ishift, jshift)
   integer,                                       intent(in)    :: is, ie, js, je
   integer,   dimension(:),                       intent(in)    :: pelist
   integer(i8_kind), dimension(:,:), contiguous, target, intent(in)    :: array_seg
   integer(i8_kind), dimension(:,:), contiguous, target, intent(inout) :: gather_data
   integer,   dimension(2),                       intent(in)    :: axis_to_storage
   logical,                                       intent(in)    :: is_root_pe
   integer,   optional,                           intent(in)    :: ishift, jshift

   integer(i8_kind), pointer ::  arr3D(:,:,:)
   integer(i8_kind), pointer :: data3D(:,:,:)

   arr3D(1:size(array_seg,1),1:size(array_seg,2),1:1) => array_seg
   if (is_root_pe) then
     data3D(1:size(gather_data,1),1:size(gather_data,2),1:1) => gather_data
   else
     data3D => null()
   endif

   call mpp_gather(is, ie, js, je, 1, pelist, arr3D, data3D, [axis_to_storage, 3], is_root_pe, &
                   ishift, jshift)
   return

end subroutine mpp_gather_pelist_int8_gen_2d


subroutine mpp_gather_pelist_int8_3d(is, ie, js, je, nk, pelist, array_seg, gather_data, is_root_pe, &
                                 ishift, jshift)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   integer(i8_kind), dimension(is:ie,js:je,1:nk), intent(in)    :: array_seg
   integer(i8_kind), dimension(:,:,:),            intent(inout) :: gather_data
   logical,                                intent(in)    :: is_root_pe
   integer,   optional,                    intent(in)    :: ishift, jshift

   integer, dimension(3) :: axis_to_storage

   axis_to_storage = (/1, 2, 3/)

   call mpp_gather(is, ie, js, je, nk, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                   ishift, jshift)
   return

end subroutine mpp_gather_pelist_int8_3d

subroutine mpp_gather_pelist_int8_gen_3d(is, ie, js, je, nk, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                                     ishift, jshift)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   integer(i8_kind), dimension(:,:,:),            intent(in)    :: array_seg
   integer(i8_kind), dimension(:,:,:),            intent(inout) :: gather_data
   integer,   dimension(3),                intent(in)    :: axis_to_storage
   logical,                                intent(in)    :: is_root_pe
   integer,   optional,                    intent(in)    :: ishift, jshift

   integer :: root_pe, root_pe_test
   integer :: k, us, ue, vs, ve, ws, we
   integer :: i1, i2, j1, j2, ioff, joff
   integer :: base_idx, send_count, msg_start
   integer :: blocksize_u, blocksize_v, blocksize_w, blocksize
   integer, dimension(3) :: start_idx, stop_idx, storage_to_axis
   integer, dimension(:), allocatable :: gind, counts
   integer(i8_kind), dimension(:), allocatable :: rbuf

   if (.not.ANY(mpp_pe().eq.pelist(:))) return

   ! Check axis_to_storage is a permutation of 1..3
   if ( any(axis_to_storage < 1) .or. any(axis_to_storage > 3) ) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): axis_to_storage entries must be in {1,2,3}")
   if ( axis_to_storage(1) == axis_to_storage(2) .or. axis_to_storage(1) == axis_to_storage(3) &
        .or. axis_to_storage(2) == axis_to_storage(3) ) &
     call mpp_error(FATAL, "fms_io(mpp_gather_pelist): axis_to_storage must be a permutation of 1,2,3")

   if (is_root_pe) then
     root_pe = mpp_pe()
     root_pe_test = 999
     if (.not.ANY(pelist(:).eq.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): root_pe not a member of pelist")
   else
     root_pe = 0
     root_pe_test = -999
   endif
! need this check in case MPI-rank 0 is a member of the pelist
   call mpp_max(root_pe_test, pelist)
   if (root_pe_test.lt.0) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): root_pe not specified or not a member of the pelist")
! need to make sure only one root_pe has been specified
   call mpp_sum(root_pe, pelist)
   if ((is_root_pe) .and. (mpp_pe().ne.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): too many root_pes specified")

   ioff=0
   joff=0
   if (present(ishift)) ioff=ishift
   if (present(jshift)) joff=jshift

   ! Initialize storage dim to axis map
   storage_to_axis(axis_to_storage(1)) = 1
   storage_to_axis(axis_to_storage(2)) = 2
   storage_to_axis(axis_to_storage(3)) = 3

   ! gather indices into global index on root_pe
   if (is_root_pe) then
     allocate(gind(4*size(pelist)))
   else
     allocate(gind(1))
   endif
   call mpp_gather((/is, ie, js, je/), gind, pelist)

   ! Compute recv counts and allocate 1d recv buffer (rbuf)
   if (is_root_pe) then
      allocate(counts(size(pelist)))

      do k = 1, size(pelist)
         base_idx = 4*(k-1)
         i1 = gind( base_idx + 1 ) ;; i2 = gind( base_idx + 2 )
         j1 = gind( base_idx + 3 ) ;; j2 = gind( base_idx + 4 )
         counts(k) = (i2 - i1 + 1) * (j2 - j1 + 1) * nk
      enddo

      allocate(rbuf(sum(counts)))
   else
      ! Non-root: MPI ignores recv args, but they must still be valid actual arguements
      allocate(counts(1))
      counts = 0
      allocate(rbuf(1))
   endif

   send_count = (ie-is+1)*(je-js+1)*nk

   ! Get generalized stop indicies for array_seg
   stop_idx = (/ie-is+1, je-js+1, nk/)
   ue = stop_idx(storage_to_axis(1))
   ve = stop_idx(storage_to_axis(2))
   we = stop_idx(storage_to_axis(3))

   ! gather data into 1d recv buffer
   call mpp_gather(reshape(array_seg(1:ue,1:ve,1:we),[send_count]), send_count, rbuf, counts, pelist)

   ! Unpack recv buffer into return array (gather_data)
   if (is_root_pe) then
      msg_start = 1
      do k = 1, size(pelist)
         base_idx = 4*(k-1)
         i1 = gind( base_idx + 1 ) + ioff ;; i2 = gind( base_idx + 2 ) + ioff
         j1 = gind( base_idx + 3 ) + joff ;; j2 = gind( base_idx + 4 ) + joff

         ! Get generalized start/stop indicies
         start_idx = (/i1,j1,1/)
         stop_idx  = (/i2,j2,nk/)

         us = start_idx(storage_to_axis(1)) ;; ue = stop_idx(storage_to_axis(1))
         vs = start_idx(storage_to_axis(2)) ;; ve = stop_idx(storage_to_axis(2))
         ws = start_idx(storage_to_axis(3)) ;; we = stop_idx(storage_to_axis(3))

         ! Compute block sizes
         blocksize_u = ue - us + 1
         blocksize_v = ve - vs + 1
         blocksize_w = we - ws + 1
         blocksize   = blocksize_u * blocksize_v * blocksize_w

         gather_data(us:ue, vs:ve, ws:we) = reshape(rbuf(msg_start:msg_start+blocksize-1), &
                                                         [blocksize_u, blocksize_v, blocksize_w])

         msg_start = msg_start + blocksize
      enddo

      deallocate(gind)
   endif

   deallocate(rbuf, counts)

   call mpp_sync_self()

end subroutine mpp_gather_pelist_int8_gen_3d
# 437 "mpp/include/mpp_comm.inc" 2


# 456 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_gather.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_gather_real4_1d(sbuf, rbuf, pelist)
! JWD: Did not create mpp_gather_2d because have no requirement for it
! JWD: See mpp_gather_2dv below
   real(r4_kind), dimension(:),    intent(in) :: sbuf
   real(r4_kind), dimension(:), intent(inout) :: rbuf
   integer,   dimension(:),    intent(in), optional :: pelist(:)

   integer :: cnt, l, nproc, op_root, ierr
   integer, allocatable :: pelist2(:)

   if( .NOT.module_is_initialized ) call mpp_error( FATAL, 'MPP_GATHER_1D_: You must first call mpp_init.' )

!  If pelist is provided, the first position must be the operation root, w.r.t. new comm, op_root = 0
   if(PRESENT(pelist))then
      if(.not.ANY(mpp_pe().eq.pelist(:))) return
      nproc = size(pelist)
      allocate(pelist2(nproc))
      pelist2 = pelist
   else
      nproc = mpp_npes()
      allocate(pelist2(nproc))
      pelist2 = (/ (l, l=root_pe, nproc-1+root_pe) /)
   endif
   op_root = pelist2(1)

   cnt = size(sbuf(:))
   if((mpp_pe().eq.op_root).AND.(size(rbuf(:)) < cnt*nproc)) call mpp_error(FATAL, &
          "MPP_GATHER_1D_: size(rbuf) must be at least npes*size(sbuf) ")

   call mpp_gather( sbuf, rbuf, size(sbuf), op_root, pelist2, ierr )

   call mpp_sync_self()
   deallocate(pelist2)
end subroutine mpp_gather_real4_1d

subroutine mpp_gather_real4_1dv(sbuf, ssize, rbuf, rsize, pelist)
   real(r4_kind), dimension(:),    intent(in) :: sbuf
   real(r4_kind), dimension(:), intent(inout) :: rbuf
   integer,                    intent(in) :: ssize
   integer,   dimension(:),    intent(in) :: rsize
   integer,   dimension(:),    intent(in), optional :: pelist(:)

   integer :: l, nproc, op_root, ierr
   integer, dimension(:), allocatable :: displs
   integer, dimension(:), allocatable :: pelist2

!  If pelist is provided, the first position must be
!  the operation root
   if(PRESENT(pelist))then
      nproc = size(pelist)
      allocate(pelist2(nproc))
      pelist2 = pelist
   else
      nproc = mpp_npes()
      allocate(pelist2(nproc))
      pelist2 = (/ (l, l=root_pe, nproc-1+root_pe) /)
   endif
   op_root = pelist2(1)

   if(pe .eq. op_root) then
      allocate(displs(nproc))

      displs(1) = 0
      do l = 2, nproc
         displs(l) = displs(l-1) + rsize(l-1)
      enddo
   else
      allocate(displs(1))
   endif

   call mpp_gather( sbuf, ssize, rbuf, rsize, displs, op_root, pelist2, ierr )

   call mpp_sync_self()
   deallocate(pelist2)
   deallocate(displs)
end subroutine mpp_gather_real4_1dv

subroutine mpp_gather_pelist_real4_2d(is, ie, js, je, pelist, array_seg, gather_data, is_root_pe, &
                                 ishift, jshift)
   integer,                                       intent(in)    :: is, ie, js, je
   integer,   dimension(:),                       intent(in)    :: pelist
   real(r4_kind), dimension(is:ie,js:je),     target, intent(in)    :: array_seg
   real(r4_kind), dimension(:,:), contiguous, target, intent(inout) :: gather_data
   logical,                                       intent(in)    :: is_root_pe
   integer,   optional,                           intent(in)    :: ishift, jshift

   integer, dimension(2)  :: axis_to_storage

   axis_to_storage = (/1,2/)

   call mpp_gather(is, ie, js, je, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                   ishift, jshift)
   return

end subroutine mpp_gather_pelist_real4_2d

subroutine mpp_gather_pelist_real4_gen_2d(is, ie, js, je, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                                     ishift, jshift)
   integer,                                       intent(in)    :: is, ie, js, je
   integer,   dimension(:),                       intent(in)    :: pelist
   real(r4_kind), dimension(:,:), contiguous, target, intent(in)    :: array_seg
   real(r4_kind), dimension(:,:), contiguous, target, intent(inout) :: gather_data
   integer,   dimension(2),                       intent(in)    :: axis_to_storage
   logical,                                       intent(in)    :: is_root_pe
   integer,   optional,                           intent(in)    :: ishift, jshift

   real(r4_kind), pointer ::  arr3D(:,:,:)
   real(r4_kind), pointer :: data3D(:,:,:)

   arr3D(1:size(array_seg,1),1:size(array_seg,2),1:1) => array_seg
   if (is_root_pe) then
     data3D(1:size(gather_data,1),1:size(gather_data,2),1:1) => gather_data
   else
     data3D => null()
   endif

   call mpp_gather(is, ie, js, je, 1, pelist, arr3D, data3D, [axis_to_storage, 3], is_root_pe, &
                   ishift, jshift)
   return

end subroutine mpp_gather_pelist_real4_gen_2d


subroutine mpp_gather_pelist_real4_3d(is, ie, js, je, nk, pelist, array_seg, gather_data, is_root_pe, &
                                 ishift, jshift)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   real(r4_kind), dimension(is:ie,js:je,1:nk), intent(in)    :: array_seg
   real(r4_kind), dimension(:,:,:),            intent(inout) :: gather_data
   logical,                                intent(in)    :: is_root_pe
   integer,   optional,                    intent(in)    :: ishift, jshift

   integer, dimension(3) :: axis_to_storage

   axis_to_storage = (/1, 2, 3/)

   call mpp_gather(is, ie, js, je, nk, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                   ishift, jshift)
   return

end subroutine mpp_gather_pelist_real4_3d

subroutine mpp_gather_pelist_real4_gen_3d(is, ie, js, je, nk, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                                     ishift, jshift)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   real(r4_kind), dimension(:,:,:),            intent(in)    :: array_seg
   real(r4_kind), dimension(:,:,:),            intent(inout) :: gather_data
   integer,   dimension(3),                intent(in)    :: axis_to_storage
   logical,                                intent(in)    :: is_root_pe
   integer,   optional,                    intent(in)    :: ishift, jshift

   integer :: root_pe, root_pe_test
   integer :: k, us, ue, vs, ve, ws, we
   integer :: i1, i2, j1, j2, ioff, joff
   integer :: base_idx, send_count, msg_start
   integer :: blocksize_u, blocksize_v, blocksize_w, blocksize
   integer, dimension(3) :: start_idx, stop_idx, storage_to_axis
   integer, dimension(:), allocatable :: gind, counts
   real(r4_kind), dimension(:), allocatable :: rbuf

   if (.not.ANY(mpp_pe().eq.pelist(:))) return

   ! Check axis_to_storage is a permutation of 1..3
   if ( any(axis_to_storage < 1) .or. any(axis_to_storage > 3) ) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): axis_to_storage entries must be in {1,2,3}")
   if ( axis_to_storage(1) == axis_to_storage(2) .or. axis_to_storage(1) == axis_to_storage(3) &
        .or. axis_to_storage(2) == axis_to_storage(3) ) &
     call mpp_error(FATAL, "fms_io(mpp_gather_pelist): axis_to_storage must be a permutation of 1,2,3")

   if (is_root_pe) then
     root_pe = mpp_pe()
     root_pe_test = 999
     if (.not.ANY(pelist(:).eq.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): root_pe not a member of pelist")
   else
     root_pe = 0
     root_pe_test = -999
   endif
! need this check in case MPI-rank 0 is a member of the pelist
   call mpp_max(root_pe_test, pelist)
   if (root_pe_test.lt.0) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): root_pe not specified or not a member of the pelist")
! need to make sure only one root_pe has been specified
   call mpp_sum(root_pe, pelist)
   if ((is_root_pe) .and. (mpp_pe().ne.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): too many root_pes specified")

   ioff=0
   joff=0
   if (present(ishift)) ioff=ishift
   if (present(jshift)) joff=jshift

   ! Initialize storage dim to axis map
   storage_to_axis(axis_to_storage(1)) = 1
   storage_to_axis(axis_to_storage(2)) = 2
   storage_to_axis(axis_to_storage(3)) = 3

   ! gather indices into global index on root_pe
   if (is_root_pe) then
     allocate(gind(4*size(pelist)))
   else
     allocate(gind(1))
   endif
   call mpp_gather((/is, ie, js, je/), gind, pelist)

   ! Compute recv counts and allocate 1d recv buffer (rbuf)
   if (is_root_pe) then
      allocate(counts(size(pelist)))

      do k = 1, size(pelist)
         base_idx = 4*(k-1)
         i1 = gind( base_idx + 1 ) ;; i2 = gind( base_idx + 2 )
         j1 = gind( base_idx + 3 ) ;; j2 = gind( base_idx + 4 )
         counts(k) = (i2 - i1 + 1) * (j2 - j1 + 1) * nk
      enddo

      allocate(rbuf(sum(counts)))
   else
      ! Non-root: MPI ignores recv args, but they must still be valid actual arguements
      allocate(counts(1))
      counts = 0
      allocate(rbuf(1))
   endif

   send_count = (ie-is+1)*(je-js+1)*nk

   ! Get generalized stop indicies for array_seg
   stop_idx = (/ie-is+1, je-js+1, nk/)
   ue = stop_idx(storage_to_axis(1))
   ve = stop_idx(storage_to_axis(2))
   we = stop_idx(storage_to_axis(3))

   ! gather data into 1d recv buffer
   call mpp_gather(reshape(array_seg(1:ue,1:ve,1:we),[send_count]), send_count, rbuf, counts, pelist)

   ! Unpack recv buffer into return array (gather_data)
   if (is_root_pe) then
      msg_start = 1
      do k = 1, size(pelist)
         base_idx = 4*(k-1)
         i1 = gind( base_idx + 1 ) + ioff ;; i2 = gind( base_idx + 2 ) + ioff
         j1 = gind( base_idx + 3 ) + joff ;; j2 = gind( base_idx + 4 ) + joff

         ! Get generalized start/stop indicies
         start_idx = (/i1,j1,1/)
         stop_idx  = (/i2,j2,nk/)

         us = start_idx(storage_to_axis(1)) ;; ue = stop_idx(storage_to_axis(1))
         vs = start_idx(storage_to_axis(2)) ;; ve = stop_idx(storage_to_axis(2))
         ws = start_idx(storage_to_axis(3)) ;; we = stop_idx(storage_to_axis(3))

         ! Compute block sizes
         blocksize_u = ue - us + 1
         blocksize_v = ve - vs + 1
         blocksize_w = we - ws + 1
         blocksize   = blocksize_u * blocksize_v * blocksize_w

         gather_data(us:ue, vs:ve, ws:we) = reshape(rbuf(msg_start:msg_start+blocksize-1), &
                                                         [blocksize_u, blocksize_v, blocksize_w])

         msg_start = msg_start + blocksize
      enddo

      deallocate(gind)
   endif

   deallocate(rbuf, counts)

   call mpp_sync_self()

end subroutine mpp_gather_pelist_real4_gen_3d
# 456 "mpp/include/mpp_comm.inc" 2

# 474 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_gather.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
subroutine mpp_gather_real8_1d(sbuf, rbuf, pelist)
! JWD: Did not create mpp_gather_2d because have no requirement for it
! JWD: See mpp_gather_2dv below
   real(r8_kind), dimension(:),    intent(in) :: sbuf
   real(r8_kind), dimension(:), intent(inout) :: rbuf
   integer,   dimension(:),    intent(in), optional :: pelist(:)

   integer :: cnt, l, nproc, op_root, ierr
   integer, allocatable :: pelist2(:)

   if( .NOT.module_is_initialized ) call mpp_error( FATAL, 'MPP_GATHER_1D_: You must first call mpp_init.' )

!  If pelist is provided, the first position must be the operation root, w.r.t. new comm, op_root = 0
   if(PRESENT(pelist))then
      if(.not.ANY(mpp_pe().eq.pelist(:))) return
      nproc = size(pelist)
      allocate(pelist2(nproc))
      pelist2 = pelist
   else
      nproc = mpp_npes()
      allocate(pelist2(nproc))
      pelist2 = (/ (l, l=root_pe, nproc-1+root_pe) /)
   endif
   op_root = pelist2(1)

   cnt = size(sbuf(:))
   if((mpp_pe().eq.op_root).AND.(size(rbuf(:)) < cnt*nproc)) call mpp_error(FATAL, &
          "MPP_GATHER_1D_: size(rbuf) must be at least npes*size(sbuf) ")

   call mpp_gather( sbuf, rbuf, size(sbuf), op_root, pelist2, ierr )

   call mpp_sync_self()
   deallocate(pelist2)
end subroutine mpp_gather_real8_1d

subroutine mpp_gather_real8_1dv(sbuf, ssize, rbuf, rsize, pelist)
   real(r8_kind), dimension(:),    intent(in) :: sbuf
   real(r8_kind), dimension(:), intent(inout) :: rbuf
   integer,                    intent(in) :: ssize
   integer,   dimension(:),    intent(in) :: rsize
   integer,   dimension(:),    intent(in), optional :: pelist(:)

   integer :: l, nproc, op_root, ierr
   integer, dimension(:), allocatable :: displs
   integer, dimension(:), allocatable :: pelist2

!  If pelist is provided, the first position must be
!  the operation root
   if(PRESENT(pelist))then
      nproc = size(pelist)
      allocate(pelist2(nproc))
      pelist2 = pelist
   else
      nproc = mpp_npes()
      allocate(pelist2(nproc))
      pelist2 = (/ (l, l=root_pe, nproc-1+root_pe) /)
   endif
   op_root = pelist2(1)

   if(pe .eq. op_root) then
      allocate(displs(nproc))

      displs(1) = 0
      do l = 2, nproc
         displs(l) = displs(l-1) + rsize(l-1)
      enddo
   else
      allocate(displs(1))
   endif

   call mpp_gather( sbuf, ssize, rbuf, rsize, displs, op_root, pelist2, ierr )

   call mpp_sync_self()
   deallocate(pelist2)
   deallocate(displs)
end subroutine mpp_gather_real8_1dv

subroutine mpp_gather_pelist_real8_2d(is, ie, js, je, pelist, array_seg, gather_data, is_root_pe, &
                                 ishift, jshift)
   integer,                                       intent(in)    :: is, ie, js, je
   integer,   dimension(:),                       intent(in)    :: pelist
   real(r8_kind), dimension(is:ie,js:je),     target, intent(in)    :: array_seg
   real(r8_kind), dimension(:,:), contiguous, target, intent(inout) :: gather_data
   logical,                                       intent(in)    :: is_root_pe
   integer,   optional,                           intent(in)    :: ishift, jshift

   integer, dimension(2)  :: axis_to_storage

   axis_to_storage = (/1,2/)

   call mpp_gather(is, ie, js, je, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                   ishift, jshift)
   return

end subroutine mpp_gather_pelist_real8_2d

subroutine mpp_gather_pelist_real8_gen_2d(is, ie, js, je, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                                     ishift, jshift)
   integer,                                       intent(in)    :: is, ie, js, je
   integer,   dimension(:),                       intent(in)    :: pelist
   real(r8_kind), dimension(:,:), contiguous, target, intent(in)    :: array_seg
   real(r8_kind), dimension(:,:), contiguous, target, intent(inout) :: gather_data
   integer,   dimension(2),                       intent(in)    :: axis_to_storage
   logical,                                       intent(in)    :: is_root_pe
   integer,   optional,                           intent(in)    :: ishift, jshift

   real(r8_kind), pointer ::  arr3D(:,:,:)
   real(r8_kind), pointer :: data3D(:,:,:)

   arr3D(1:size(array_seg,1),1:size(array_seg,2),1:1) => array_seg
   if (is_root_pe) then
     data3D(1:size(gather_data,1),1:size(gather_data,2),1:1) => gather_data
   else
     data3D => null()
   endif

   call mpp_gather(is, ie, js, je, 1, pelist, arr3D, data3D, [axis_to_storage, 3], is_root_pe, &
                   ishift, jshift)
   return

end subroutine mpp_gather_pelist_real8_gen_2d


subroutine mpp_gather_pelist_real8_3d(is, ie, js, je, nk, pelist, array_seg, gather_data, is_root_pe, &
                                 ishift, jshift)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   real(r8_kind), dimension(is:ie,js:je,1:nk), intent(in)    :: array_seg
   real(r8_kind), dimension(:,:,:),            intent(inout) :: gather_data
   logical,                                intent(in)    :: is_root_pe
   integer,   optional,                    intent(in)    :: ishift, jshift

   integer, dimension(3) :: axis_to_storage

   axis_to_storage = (/1, 2, 3/)

   call mpp_gather(is, ie, js, je, nk, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                   ishift, jshift)
   return

end subroutine mpp_gather_pelist_real8_3d

subroutine mpp_gather_pelist_real8_gen_3d(is, ie, js, je, nk, pelist, array_seg, gather_data, axis_to_storage, is_root_pe, &
                                     ishift, jshift)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   real(r8_kind), dimension(:,:,:),            intent(in)    :: array_seg
   real(r8_kind), dimension(:,:,:),            intent(inout) :: gather_data
   integer,   dimension(3),                intent(in)    :: axis_to_storage
   logical,                                intent(in)    :: is_root_pe
   integer,   optional,                    intent(in)    :: ishift, jshift

   integer :: root_pe, root_pe_test
   integer :: k, us, ue, vs, ve, ws, we
   integer :: i1, i2, j1, j2, ioff, joff
   integer :: base_idx, send_count, msg_start
   integer :: blocksize_u, blocksize_v, blocksize_w, blocksize
   integer, dimension(3) :: start_idx, stop_idx, storage_to_axis
   integer, dimension(:), allocatable :: gind, counts
   real(r8_kind), dimension(:), allocatable :: rbuf

   if (.not.ANY(mpp_pe().eq.pelist(:))) return

   ! Check axis_to_storage is a permutation of 1..3
   if ( any(axis_to_storage < 1) .or. any(axis_to_storage > 3) ) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): axis_to_storage entries must be in {1,2,3}")
   if ( axis_to_storage(1) == axis_to_storage(2) .or. axis_to_storage(1) == axis_to_storage(3) &
        .or. axis_to_storage(2) == axis_to_storage(3) ) &
     call mpp_error(FATAL, "fms_io(mpp_gather_pelist): axis_to_storage must be a permutation of 1,2,3")

   if (is_root_pe) then
     root_pe = mpp_pe()
     root_pe_test = 999
     if (.not.ANY(pelist(:).eq.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): root_pe not a member of pelist")
   else
     root_pe = 0
     root_pe_test = -999
   endif
! need this check in case MPI-rank 0 is a member of the pelist
   call mpp_max(root_pe_test, pelist)
   if (root_pe_test.lt.0) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): root_pe not specified or not a member of the pelist")
! need to make sure only one root_pe has been specified
   call mpp_sum(root_pe, pelist)
   if ((is_root_pe) .and. (mpp_pe().ne.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_gather_pelist): too many root_pes specified")

   ioff=0
   joff=0
   if (present(ishift)) ioff=ishift
   if (present(jshift)) joff=jshift

   ! Initialize storage dim to axis map
   storage_to_axis(axis_to_storage(1)) = 1
   storage_to_axis(axis_to_storage(2)) = 2
   storage_to_axis(axis_to_storage(3)) = 3

   ! gather indices into global index on root_pe
   if (is_root_pe) then
     allocate(gind(4*size(pelist)))
   else
     allocate(gind(1))
   endif
   call mpp_gather((/is, ie, js, je/), gind, pelist)

   ! Compute recv counts and allocate 1d recv buffer (rbuf)
   if (is_root_pe) then
      allocate(counts(size(pelist)))

      do k = 1, size(pelist)
         base_idx = 4*(k-1)
         i1 = gind( base_idx + 1 ) ;; i2 = gind( base_idx + 2 )
         j1 = gind( base_idx + 3 ) ;; j2 = gind( base_idx + 4 )
         counts(k) = (i2 - i1 + 1) * (j2 - j1 + 1) * nk
      enddo

      allocate(rbuf(sum(counts)))
   else
      ! Non-root: MPI ignores recv args, but they must still be valid actual arguements
      allocate(counts(1))
      counts = 0
      allocate(rbuf(1))
   endif

   send_count = (ie-is+1)*(je-js+1)*nk

   ! Get generalized stop indicies for array_seg
   stop_idx = (/ie-is+1, je-js+1, nk/)
   ue = stop_idx(storage_to_axis(1))
   ve = stop_idx(storage_to_axis(2))
   we = stop_idx(storage_to_axis(3))

   ! gather data into 1d recv buffer
   call mpp_gather(reshape(array_seg(1:ue,1:ve,1:we),[send_count]), send_count, rbuf, counts, pelist)

   ! Unpack recv buffer into return array (gather_data)
   if (is_root_pe) then
      msg_start = 1
      do k = 1, size(pelist)
         base_idx = 4*(k-1)
         i1 = gind( base_idx + 1 ) + ioff ;; i2 = gind( base_idx + 2 ) + ioff
         j1 = gind( base_idx + 3 ) + joff ;; j2 = gind( base_idx + 4 ) + joff

         ! Get generalized start/stop indicies
         start_idx = (/i1,j1,1/)
         stop_idx  = (/i2,j2,nk/)

         us = start_idx(storage_to_axis(1)) ;; ue = stop_idx(storage_to_axis(1))
         vs = start_idx(storage_to_axis(2)) ;; ve = stop_idx(storage_to_axis(2))
         ws = start_idx(storage_to_axis(3)) ;; we = stop_idx(storage_to_axis(3))

         ! Compute block sizes
         blocksize_u = ue - us + 1
         blocksize_v = ve - vs + 1
         blocksize_w = we - ws + 1
         blocksize   = blocksize_u * blocksize_v * blocksize_w

         gather_data(us:ue, vs:ve, ws:we) = reshape(rbuf(msg_start:msg_start+blocksize-1), &
                                                         [blocksize_u, blocksize_v, blocksize_w])

         msg_start = msg_start + blocksize
      enddo

      deallocate(gind)
   endif

   deallocate(rbuf, counts)

   call mpp_sync_self()

end subroutine mpp_gather_pelist_real8_gen_3d
# 474 "mpp/include/mpp_comm.inc" 2

!#################################################
# 489 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_scatter.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> addtogroup mpp_mod

!> @brief Scatter data from one pe to the specified pes.
!!
!> Scatter (ie - is) * (je - js) contiguous elements of array data from the designated root pe
!! into contigous members of array segment in each pe that is included in the pelist argument.
subroutine mpp_scatter_pelist_int4_2d(is, ie, js, je, pelist, array_seg, input_data, is_root_pe)
   integer,                           intent(in)    :: is, ie, js, je !< indices of segment array
   integer,   dimension(:),           intent(in)    :: pelist !<PE list of target pes,
                                                              !! must be in monotonic increasing order
   integer(i4_kind), dimension(is:ie,js:je), target, intent(inout)  :: array_seg !< 2D array of output data
   integer(i4_kind), dimension(:,:), contiguous, target, intent(in) :: input_data !< 2D array of input data
   logical,                           intent(in)    :: is_root_pe !< operational root pe

   integer, dimension(2)   :: axis_to_storage

   axis_to_storage = (/1,2/)

   call mpp_scatter(is, ie, js, je, pelist, array_seg, input_data, axis_to_storage, is_root_pe)

   return

end subroutine mpp_scatter_pelist_int4_2d

subroutine mpp_scatter_pelist_int4_gen_2d(is, ie, js, je, pelist, array_seg, input_data, axis_to_storage, is_root_pe)
   integer,                           intent(in)    :: is, ie, js, je !< indices of segment array
   integer,   dimension(:),           intent(in)    :: pelist !<PE list of target pes,
                                                              !! must be in monotonic increasing order
   integer(i4_kind), dimension(:,:), contiguous, target, intent(inout)  :: array_seg !< 2D array of output data
   integer(i4_kind), dimension(:,:), contiguous, target, intent(in) :: input_data !< 2D array of input data
   integer, dimension(2),             intent(in)    :: axis_to_storage
   logical,                           intent(in)    :: is_root_pe !< operational root pe

   integer(i4_kind), pointer ::  arr3D(:,:,:)
   integer(i4_kind), pointer :: data3D(:,:,:)

   arr3D(1:size(array_seg,1),1:size(array_seg,2),1:1) => array_seg
   if (is_root_pe) then
     data3D(1:size(input_data,1),1:size(input_data,2),1:1) => input_data
   else
     data3D => null()
   endif

   call mpp_scatter(is, ie, js, je, 1, pelist, arr3D, data3D, [axis_to_storage, 3], is_root_pe)

   return

end subroutine mpp_scatter_pelist_int4_gen_2d

subroutine mpp_scatter_pelist_int4_3d(is, ie, js, je, nk, pelist, array_seg, input_data, is_root_pe)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   integer(i4_kind), dimension(:,:,:),            intent(inout) :: array_seg
   integer(i4_kind), dimension(:,:,:),            intent(in)    :: input_data
   logical,                                intent(in)    :: is_root_pe

   integer, dimension(3) :: axis_to_storage

   axis_to_storage = (/1,2,3/)

   call mpp_scatter(is, ie, js, je, nk, pelist, array_seg, input_data, axis_to_storage, is_root_pe)

   return

end subroutine mpp_scatter_pelist_int4_3d

subroutine mpp_scatter_pelist_int4_gen_3d(is, ie, js, je, nk, pelist, array_seg, input_data, axis_to_storage, is_root_pe)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   integer(i4_kind), dimension(:,:,:),            intent(inout) :: array_seg
   integer(i4_kind), dimension(:,:,:),            intent(in)    :: input_data
   integer, dimension(3),                  intent(in)    :: axis_to_storage
   logical,                                intent(in)    :: is_root_pe

   integer :: i, j, k, n, m, ierr, base_idx
   integer :: i1, i2, j1, j2
   integer :: us, ue, vs, ve, ws, we
   integer :: root_pe, root_pe_test, recv_count
   integer, dimension(size(pelist)) :: counts, displs
   integer, dimension(4*size(pelist)) :: gind
   integer, dimension(3) :: start_idx, end_idx, storage_to_axis
   integer(i4_kind), dimension(:), allocatable :: temp

   if (.not.ANY(mpp_pe().eq.pelist(:))) return

   if (is_root_pe) then
     root_pe = mpp_pe()
     root_pe_test = 999
     if (.not.ANY(pelist(:).eq.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_scatter_pelist): root_pe not a member of pelist")
   else
     root_pe = 0
     root_pe_test = -999
   endif

! need this check in case MPI-rank 0 is a member of the pelist
   call mpp_max(root_pe_test, pelist)
   if (root_pe_test.lt.0) call mpp_error(FATAL, &
                "fms_io(mpp_scatter_pelist): root_pe not specified or not a member of the pelist")
! need to make sure only one root_pe has been specified
   call mpp_sum(root_pe, pelist)
   if ((is_root_pe) .and. (mpp_pe().ne.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_scatter_pelist): too many root_pes specified")

   ! Initialize storage dim to axis map
   storage_to_axis(axis_to_storage(1)) = 1
   storage_to_axis(axis_to_storage(2)) = 2
   storage_to_axis(axis_to_storage(3)) = 3


   ! Gather Indices on root pe
   call mpp_gather((/is, ie, js, je/), gind, pelist)

   ! Compute counts, displs, and setup 1d send buffer (temp)
   if (is_root_pe) then

      displs(1) = 0
      do k = 1, size(pelist)
         base_idx = 4*(k-1)
         i1 = gind( base_idx + 1 ) ;; i2 = gind( base_idx + 2 )
         j1 = gind( base_idx + 3 ) ;; j2 = gind( base_idx + 4 )
         counts(k) = (i2 - i1 + 1) * (j2 - j1 + 1) * nk
         if (k > 1) displs(k) = displs(k-1) + counts(k-1)
      enddo

      allocate(temp(sum(counts)))

      m = 1
      do n = 1, size(pelist)
         base_idx = 4*(n-1)
         i1 = gind( base_idx + 1 ) ;; i2 = gind( base_idx + 2 )
         j1 = gind( base_idx + 3 ) ;; j2 = gind( base_idx + 4 )

         start_idx = (/i1, j1, 1/)
         end_idx   = (/i2, j2, nk/)

         us = start_idx(storage_to_axis(1))  ;;  ue = end_idx(storage_to_axis(1))
         vs = start_idx(storage_to_axis(2))  ;;  ve = end_idx(storage_to_axis(2))
         ws = start_idx(storage_to_axis(3))  ;;  we = end_idx(storage_to_axis(3))

         temp(m:m+counts(n)-1) = reshape( input_data(us:ue, vs:ve, ws:we), [counts(n)] )
         m = m + counts(n)
       enddo
   else
      allocate(temp(1))
   endif

   ! Compute recv_count on each rank
   recv_count = (ie-is+1)*(je-js+1)*nk

   call mpp_scatter(temp, counts, displs, array_seg, recv_count, root_pe, pelist, ierr)

   call mpp_sync_self()

   deallocate(temp)

   return

end subroutine mpp_scatter_pelist_int4_gen_3d
# 489 "mpp/include/mpp_comm.inc" 2

# 503 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_scatter.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> addtogroup mpp_mod

!> @brief Scatter data from one pe to the specified pes.
!!
!> Scatter (ie - is) * (je - js) contiguous elements of array data from the designated root pe
!! into contigous members of array segment in each pe that is included in the pelist argument.
subroutine mpp_scatter_pelist_int8_2d(is, ie, js, je, pelist, array_seg, input_data, is_root_pe)
   integer,                           intent(in)    :: is, ie, js, je !< indices of segment array
   integer,   dimension(:),           intent(in)    :: pelist !<PE list of target pes,
                                                              !! must be in monotonic increasing order
   integer(i8_kind), dimension(is:ie,js:je), target, intent(inout)  :: array_seg !< 2D array of output data
   integer(i8_kind), dimension(:,:), contiguous, target, intent(in) :: input_data !< 2D array of input data
   logical,                           intent(in)    :: is_root_pe !< operational root pe

   integer, dimension(2)   :: axis_to_storage

   axis_to_storage = (/1,2/)

   call mpp_scatter(is, ie, js, je, pelist, array_seg, input_data, axis_to_storage, is_root_pe)

   return

end subroutine mpp_scatter_pelist_int8_2d

subroutine mpp_scatter_pelist_int8_gen_2d(is, ie, js, je, pelist, array_seg, input_data, axis_to_storage, is_root_pe)
   integer,                           intent(in)    :: is, ie, js, je !< indices of segment array
   integer,   dimension(:),           intent(in)    :: pelist !<PE list of target pes,
                                                              !! must be in monotonic increasing order
   integer(i8_kind), dimension(:,:), contiguous, target, intent(inout)  :: array_seg !< 2D array of output data
   integer(i8_kind), dimension(:,:), contiguous, target, intent(in) :: input_data !< 2D array of input data
   integer, dimension(2),             intent(in)    :: axis_to_storage
   logical,                           intent(in)    :: is_root_pe !< operational root pe

   integer(i8_kind), pointer ::  arr3D(:,:,:)
   integer(i8_kind), pointer :: data3D(:,:,:)

   arr3D(1:size(array_seg,1),1:size(array_seg,2),1:1) => array_seg
   if (is_root_pe) then
     data3D(1:size(input_data,1),1:size(input_data,2),1:1) => input_data
   else
     data3D => null()
   endif

   call mpp_scatter(is, ie, js, je, 1, pelist, arr3D, data3D, [axis_to_storage, 3], is_root_pe)

   return

end subroutine mpp_scatter_pelist_int8_gen_2d

subroutine mpp_scatter_pelist_int8_3d(is, ie, js, je, nk, pelist, array_seg, input_data, is_root_pe)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   integer(i8_kind), dimension(:,:,:),            intent(inout) :: array_seg
   integer(i8_kind), dimension(:,:,:),            intent(in)    :: input_data
   logical,                                intent(in)    :: is_root_pe

   integer, dimension(3) :: axis_to_storage

   axis_to_storage = (/1,2,3/)

   call mpp_scatter(is, ie, js, je, nk, pelist, array_seg, input_data, axis_to_storage, is_root_pe)

   return

end subroutine mpp_scatter_pelist_int8_3d

subroutine mpp_scatter_pelist_int8_gen_3d(is, ie, js, je, nk, pelist, array_seg, input_data, axis_to_storage, is_root_pe)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   integer(i8_kind), dimension(:,:,:),            intent(inout) :: array_seg
   integer(i8_kind), dimension(:,:,:),            intent(in)    :: input_data
   integer, dimension(3),                  intent(in)    :: axis_to_storage
   logical,                                intent(in)    :: is_root_pe

   integer :: i, j, k, n, m, ierr, base_idx
   integer :: i1, i2, j1, j2
   integer :: us, ue, vs, ve, ws, we
   integer :: root_pe, root_pe_test, recv_count
   integer, dimension(size(pelist)) :: counts, displs
   integer, dimension(4*size(pelist)) :: gind
   integer, dimension(3) :: start_idx, end_idx, storage_to_axis
   integer(i8_kind), dimension(:), allocatable :: temp

   if (.not.ANY(mpp_pe().eq.pelist(:))) return

   if (is_root_pe) then
     root_pe = mpp_pe()
     root_pe_test = 999
     if (.not.ANY(pelist(:).eq.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_scatter_pelist): root_pe not a member of pelist")
   else
     root_pe = 0
     root_pe_test = -999
   endif

! need this check in case MPI-rank 0 is a member of the pelist
   call mpp_max(root_pe_test, pelist)
   if (root_pe_test.lt.0) call mpp_error(FATAL, &
                "fms_io(mpp_scatter_pelist): root_pe not specified or not a member of the pelist")
! need to make sure only one root_pe has been specified
   call mpp_sum(root_pe, pelist)
   if ((is_root_pe) .and. (mpp_pe().ne.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_scatter_pelist): too many root_pes specified")

   ! Initialize storage dim to axis map
   storage_to_axis(axis_to_storage(1)) = 1
   storage_to_axis(axis_to_storage(2)) = 2
   storage_to_axis(axis_to_storage(3)) = 3


   ! Gather Indices on root pe
   call mpp_gather((/is, ie, js, je/), gind, pelist)

   ! Compute counts, displs, and setup 1d send buffer (temp)
   if (is_root_pe) then

      displs(1) = 0
      do k = 1, size(pelist)
         base_idx = 4*(k-1)
         i1 = gind( base_idx + 1 ) ;; i2 = gind( base_idx + 2 )
         j1 = gind( base_idx + 3 ) ;; j2 = gind( base_idx + 4 )
         counts(k) = (i2 - i1 + 1) * (j2 - j1 + 1) * nk
         if (k > 1) displs(k) = displs(k-1) + counts(k-1)
      enddo

      allocate(temp(sum(counts)))

      m = 1
      do n = 1, size(pelist)
         base_idx = 4*(n-1)
         i1 = gind( base_idx + 1 ) ;; i2 = gind( base_idx + 2 )
         j1 = gind( base_idx + 3 ) ;; j2 = gind( base_idx + 4 )

         start_idx = (/i1, j1, 1/)
         end_idx   = (/i2, j2, nk/)

         us = start_idx(storage_to_axis(1))  ;;  ue = end_idx(storage_to_axis(1))
         vs = start_idx(storage_to_axis(2))  ;;  ve = end_idx(storage_to_axis(2))
         ws = start_idx(storage_to_axis(3))  ;;  we = end_idx(storage_to_axis(3))

         temp(m:m+counts(n)-1) = reshape( input_data(us:ue, vs:ve, ws:we), [counts(n)] )
         m = m + counts(n)
       enddo
   else
      allocate(temp(1))
   endif

   ! Compute recv_count on each rank
   recv_count = (ie-is+1)*(je-js+1)*nk

   call mpp_scatter(temp, counts, displs, array_seg, recv_count, root_pe, pelist, ierr)

   call mpp_sync_self()

   deallocate(temp)

   return

end subroutine mpp_scatter_pelist_int8_gen_3d
# 503 "mpp/include/mpp_comm.inc" 2

# 517 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_scatter.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> addtogroup mpp_mod

!> @brief Scatter data from one pe to the specified pes.
!!
!> Scatter (ie - is) * (je - js) contiguous elements of array data from the designated root pe
!! into contigous members of array segment in each pe that is included in the pelist argument.
subroutine mpp_scatter_pelist_real4_2d(is, ie, js, je, pelist, array_seg, input_data, is_root_pe)
   integer,                           intent(in)    :: is, ie, js, je !< indices of segment array
   integer,   dimension(:),           intent(in)    :: pelist !<PE list of target pes,
                                                              !! must be in monotonic increasing order
   real(r4_kind), dimension(is:ie,js:je), target, intent(inout)  :: array_seg !< 2D array of output data
   real(r4_kind), dimension(:,:), contiguous, target, intent(in) :: input_data !< 2D array of input data
   logical,                           intent(in)    :: is_root_pe !< operational root pe

   integer, dimension(2)   :: axis_to_storage

   axis_to_storage = (/1,2/)

   call mpp_scatter(is, ie, js, je, pelist, array_seg, input_data, axis_to_storage, is_root_pe)

   return

end subroutine mpp_scatter_pelist_real4_2d

subroutine mpp_scatter_pelist_real4_gen_2d(is, ie, js, je, pelist, array_seg, input_data, axis_to_storage, is_root_pe)
   integer,                           intent(in)    :: is, ie, js, je !< indices of segment array
   integer,   dimension(:),           intent(in)    :: pelist !<PE list of target pes,
                                                              !! must be in monotonic increasing order
   real(r4_kind), dimension(:,:), contiguous, target, intent(inout)  :: array_seg !< 2D array of output data
   real(r4_kind), dimension(:,:), contiguous, target, intent(in) :: input_data !< 2D array of input data
   integer, dimension(2),             intent(in)    :: axis_to_storage
   logical,                           intent(in)    :: is_root_pe !< operational root pe

   real(r4_kind), pointer ::  arr3D(:,:,:)
   real(r4_kind), pointer :: data3D(:,:,:)

   arr3D(1:size(array_seg,1),1:size(array_seg,2),1:1) => array_seg
   if (is_root_pe) then
     data3D(1:size(input_data,1),1:size(input_data,2),1:1) => input_data
   else
     data3D => null()
   endif

   call mpp_scatter(is, ie, js, je, 1, pelist, arr3D, data3D, [axis_to_storage, 3], is_root_pe)

   return

end subroutine mpp_scatter_pelist_real4_gen_2d

subroutine mpp_scatter_pelist_real4_3d(is, ie, js, je, nk, pelist, array_seg, input_data, is_root_pe)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   real(r4_kind), dimension(:,:,:),            intent(inout) :: array_seg
   real(r4_kind), dimension(:,:,:),            intent(in)    :: input_data
   logical,                                intent(in)    :: is_root_pe

   integer, dimension(3) :: axis_to_storage

   axis_to_storage = (/1,2,3/)

   call mpp_scatter(is, ie, js, je, nk, pelist, array_seg, input_data, axis_to_storage, is_root_pe)

   return

end subroutine mpp_scatter_pelist_real4_3d

subroutine mpp_scatter_pelist_real4_gen_3d(is, ie, js, je, nk, pelist, array_seg, input_data, axis_to_storage, is_root_pe)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   real(r4_kind), dimension(:,:,:),            intent(inout) :: array_seg
   real(r4_kind), dimension(:,:,:),            intent(in)    :: input_data
   integer, dimension(3),                  intent(in)    :: axis_to_storage
   logical,                                intent(in)    :: is_root_pe

   integer :: i, j, k, n, m, ierr, base_idx
   integer :: i1, i2, j1, j2
   integer :: us, ue, vs, ve, ws, we
   integer :: root_pe, root_pe_test, recv_count
   integer, dimension(size(pelist)) :: counts, displs
   integer, dimension(4*size(pelist)) :: gind
   integer, dimension(3) :: start_idx, end_idx, storage_to_axis
   real(r4_kind), dimension(:), allocatable :: temp

   if (.not.ANY(mpp_pe().eq.pelist(:))) return

   if (is_root_pe) then
     root_pe = mpp_pe()
     root_pe_test = 999
     if (.not.ANY(pelist(:).eq.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_scatter_pelist): root_pe not a member of pelist")
   else
     root_pe = 0
     root_pe_test = -999
   endif

! need this check in case MPI-rank 0 is a member of the pelist
   call mpp_max(root_pe_test, pelist)
   if (root_pe_test.lt.0) call mpp_error(FATAL, &
                "fms_io(mpp_scatter_pelist): root_pe not specified or not a member of the pelist")
! need to make sure only one root_pe has been specified
   call mpp_sum(root_pe, pelist)
   if ((is_root_pe) .and. (mpp_pe().ne.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_scatter_pelist): too many root_pes specified")

   ! Initialize storage dim to axis map
   storage_to_axis(axis_to_storage(1)) = 1
   storage_to_axis(axis_to_storage(2)) = 2
   storage_to_axis(axis_to_storage(3)) = 3


   ! Gather Indices on root pe
   call mpp_gather((/is, ie, js, je/), gind, pelist)

   ! Compute counts, displs, and setup 1d send buffer (temp)
   if (is_root_pe) then

      displs(1) = 0
      do k = 1, size(pelist)
         base_idx = 4*(k-1)
         i1 = gind( base_idx + 1 ) ;; i2 = gind( base_idx + 2 )
         j1 = gind( base_idx + 3 ) ;; j2 = gind( base_idx + 4 )
         counts(k) = (i2 - i1 + 1) * (j2 - j1 + 1) * nk
         if (k > 1) displs(k) = displs(k-1) + counts(k-1)
      enddo

      allocate(temp(sum(counts)))

      m = 1
      do n = 1, size(pelist)
         base_idx = 4*(n-1)
         i1 = gind( base_idx + 1 ) ;; i2 = gind( base_idx + 2 )
         j1 = gind( base_idx + 3 ) ;; j2 = gind( base_idx + 4 )

         start_idx = (/i1, j1, 1/)
         end_idx   = (/i2, j2, nk/)

         us = start_idx(storage_to_axis(1))  ;;  ue = end_idx(storage_to_axis(1))
         vs = start_idx(storage_to_axis(2))  ;;  ve = end_idx(storage_to_axis(2))
         ws = start_idx(storage_to_axis(3))  ;;  we = end_idx(storage_to_axis(3))

         temp(m:m+counts(n)-1) = reshape( input_data(us:ue, vs:ve, ws:we), [counts(n)] )
         m = m + counts(n)
       enddo
   else
      allocate(temp(1))
   endif

   ! Compute recv_count on each rank
   recv_count = (ie-is+1)*(je-js+1)*nk

   call mpp_scatter(temp, counts, displs, array_seg, recv_count, root_pe, pelist, ierr)

   call mpp_sync_self()

   deallocate(temp)

   return

end subroutine mpp_scatter_pelist_real4_gen_3d
# 517 "mpp/include/mpp_comm.inc" 2

# 531 "mpp/include/mpp_comm.inc"
# 1 "mpp/include/mpp_scatter.fh" 1
!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> addtogroup mpp_mod

!> @brief Scatter data from one pe to the specified pes.
!!
!> Scatter (ie - is) * (je - js) contiguous elements of array data from the designated root pe
!! into contigous members of array segment in each pe that is included in the pelist argument.
subroutine mpp_scatter_pelist_real8_2d(is, ie, js, je, pelist, array_seg, input_data, is_root_pe)
   integer,                           intent(in)    :: is, ie, js, je !< indices of segment array
   integer,   dimension(:),           intent(in)    :: pelist !<PE list of target pes,
                                                              !! must be in monotonic increasing order
   real(r8_kind), dimension(is:ie,js:je), target, intent(inout)  :: array_seg !< 2D array of output data
   real(r8_kind), dimension(:,:), contiguous, target, intent(in) :: input_data !< 2D array of input data
   logical,                           intent(in)    :: is_root_pe !< operational root pe

   integer, dimension(2)   :: axis_to_storage

   axis_to_storage = (/1,2/)

   call mpp_scatter(is, ie, js, je, pelist, array_seg, input_data, axis_to_storage, is_root_pe)

   return

end subroutine mpp_scatter_pelist_real8_2d

subroutine mpp_scatter_pelist_real8_gen_2d(is, ie, js, je, pelist, array_seg, input_data, axis_to_storage, is_root_pe)
   integer,                           intent(in)    :: is, ie, js, je !< indices of segment array
   integer,   dimension(:),           intent(in)    :: pelist !<PE list of target pes,
                                                              !! must be in monotonic increasing order
   real(r8_kind), dimension(:,:), contiguous, target, intent(inout)  :: array_seg !< 2D array of output data
   real(r8_kind), dimension(:,:), contiguous, target, intent(in) :: input_data !< 2D array of input data
   integer, dimension(2),             intent(in)    :: axis_to_storage
   logical,                           intent(in)    :: is_root_pe !< operational root pe

   real(r8_kind), pointer ::  arr3D(:,:,:)
   real(r8_kind), pointer :: data3D(:,:,:)

   arr3D(1:size(array_seg,1),1:size(array_seg,2),1:1) => array_seg
   if (is_root_pe) then
     data3D(1:size(input_data,1),1:size(input_data,2),1:1) => input_data
   else
     data3D => null()
   endif

   call mpp_scatter(is, ie, js, je, 1, pelist, arr3D, data3D, [axis_to_storage, 3], is_root_pe)

   return

end subroutine mpp_scatter_pelist_real8_gen_2d

subroutine mpp_scatter_pelist_real8_3d(is, ie, js, je, nk, pelist, array_seg, input_data, is_root_pe)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   real(r8_kind), dimension(:,:,:),            intent(inout) :: array_seg
   real(r8_kind), dimension(:,:,:),            intent(in)    :: input_data
   logical,                                intent(in)    :: is_root_pe

   integer, dimension(3) :: axis_to_storage

   axis_to_storage = (/1,2,3/)

   call mpp_scatter(is, ie, js, je, nk, pelist, array_seg, input_data, axis_to_storage, is_root_pe)

   return

end subroutine mpp_scatter_pelist_real8_3d

subroutine mpp_scatter_pelist_real8_gen_3d(is, ie, js, je, nk, pelist, array_seg, input_data, axis_to_storage, is_root_pe)
   integer,                                intent(in)    :: is, ie, js, je, nk
   integer,   dimension(:),                intent(in)    :: pelist
   real(r8_kind), dimension(:,:,:),            intent(inout) :: array_seg
   real(r8_kind), dimension(:,:,:),            intent(in)    :: input_data
   integer, dimension(3),                  intent(in)    :: axis_to_storage
   logical,                                intent(in)    :: is_root_pe

   integer :: i, j, k, n, m, ierr, base_idx
   integer :: i1, i2, j1, j2
   integer :: us, ue, vs, ve, ws, we
   integer :: root_pe, root_pe_test, recv_count
   integer, dimension(size(pelist)) :: counts, displs
   integer, dimension(4*size(pelist)) :: gind
   integer, dimension(3) :: start_idx, end_idx, storage_to_axis
   real(r8_kind), dimension(:), allocatable :: temp

   if (.not.ANY(mpp_pe().eq.pelist(:))) return

   if (is_root_pe) then
     root_pe = mpp_pe()
     root_pe_test = 999
     if (.not.ANY(pelist(:).eq.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_scatter_pelist): root_pe not a member of pelist")
   else
     root_pe = 0
     root_pe_test = -999
   endif

! need this check in case MPI-rank 0 is a member of the pelist
   call mpp_max(root_pe_test, pelist)
   if (root_pe_test.lt.0) call mpp_error(FATAL, &
                "fms_io(mpp_scatter_pelist): root_pe not specified or not a member of the pelist")
! need to make sure only one root_pe has been specified
   call mpp_sum(root_pe, pelist)
   if ((is_root_pe) .and. (mpp_pe().ne.root_pe)) call mpp_error(FATAL, &
                "fms_io(mpp_scatter_pelist): too many root_pes specified")

   ! Initialize storage dim to axis map
   storage_to_axis(axis_to_storage(1)) = 1
   storage_to_axis(axis_to_storage(2)) = 2
   storage_to_axis(axis_to_storage(3)) = 3


   ! Gather Indices on root pe
   call mpp_gather((/is, ie, js, je/), gind, pelist)

   ! Compute counts, displs, and setup 1d send buffer (temp)
   if (is_root_pe) then

      displs(1) = 0
      do k = 1, size(pelist)
         base_idx = 4*(k-1)
         i1 = gind( base_idx + 1 ) ;; i2 = gind( base_idx + 2 )
         j1 = gind( base_idx + 3 ) ;; j2 = gind( base_idx + 4 )
         counts(k) = (i2 - i1 + 1) * (j2 - j1 + 1) * nk
         if (k > 1) displs(k) = displs(k-1) + counts(k-1)
      enddo

      allocate(temp(sum(counts)))

      m = 1
      do n = 1, size(pelist)
         base_idx = 4*(n-1)
         i1 = gind( base_idx + 1 ) ;; i2 = gind( base_idx + 2 )
         j1 = gind( base_idx + 3 ) ;; j2 = gind( base_idx + 4 )

         start_idx = (/i1, j1, 1/)
         end_idx   = (/i2, j2, nk/)

         us = start_idx(storage_to_axis(1))  ;;  ue = end_idx(storage_to_axis(1))
         vs = start_idx(storage_to_axis(2))  ;;  ve = end_idx(storage_to_axis(2))
         ws = start_idx(storage_to_axis(3))  ;;  we = end_idx(storage_to_axis(3))

         temp(m:m+counts(n)-1) = reshape( input_data(us:ue, vs:ve, ws:we), [counts(n)] )
         m = m + counts(n)
       enddo
   else
      allocate(temp(1))
   endif

   ! Compute recv_count on each rank
   recv_count = (ie-is+1)*(je-js+1)*nk

   call mpp_scatter(temp, counts, displs, array_seg, recv_count, root_pe, pelist, ierr)

   call mpp_sync_self()

   deallocate(temp)

   return

end subroutine mpp_scatter_pelist_real8_gen_3d
# 531 "mpp/include/mpp_comm.inc" 2
# 1399 "mpp/mpp.F90" 2

  end module mpp_mod
!> @}
! close documentation grouping
