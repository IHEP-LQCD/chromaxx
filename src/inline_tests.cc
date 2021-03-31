#include "meas/inline/abs_inline_measurement_factory.h"
#include "inline_tests.h"
#include "meas/glue/mesplq.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "util/ft/sftmom.h"
#include <time.h>
#include "io_general_class.h"

namespace Chroma
{
    namespace InlineTestsEnv
    {
        namespace
        {
            AbsInlineMeasurement *createMeasurement(XMLReader &xml_in,
                                                    const std::string &path)
            {
                return new Inlinetests(InlinetestsParams(xml_in, path));
            }

            //! Local registration flag
            bool registered = false;
        } // namespace

        const std::string name = "tests";

        //! Register all the factories
        bool registerAll()
        {
            bool success = true;
            if (!registered)
            {
                success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
                registered = true;
            }
            return success;
        }
    } // namespace InlineTestsEnv

    // read params
    void read(XMLReader &xml, const std::string &path, const InlinetestsParams::Param_t &param)
    {
        XMLReader paramtop(xml, path);
    }

    // write params
    void write(XMLWriter &xml, const std::string &path, const InlinetestsParams::Param_t &param)
    {
        push(xml, path);
        pop(xml);
    }

    //! gauge input
    void read(XMLReader &xml, const std::string &path, InlinetestsParams::NamedObject_t &param)
    {
        XMLReader paramtop(xml, path);

        read(paramtop, "gauge_id", param.gauge_id);
    }

    //! gauge output
    void write(XMLWriter &xml, const std::string &path, const InlinetestsParams::NamedObject_t &param)
    {
        push(xml, path);

        write(xml, "gauge_id", param.gauge_id);

        pop(xml);
    }

    // Param stuff
    InlinetestsParams::InlinetestsParams()
    {
        frequency = 0;
    }

    InlinetestsParams::InlinetestsParams(XMLReader &xml_in, const std::string &path)
    {
        try
        {
            XMLReader paramtop(xml_in, path);

            read(paramtop, "Frequency", frequency);

            // Ids
            read(paramtop, "NamedObject", named_obj);

            // Possible alternate XML file pattern
            if (paramtop.count("xml_file") != 0)
            {
                read(paramtop, "xml_file", xml_file);
            }
        }
        catch (const std::string &e)
        {
            QDPIO::cerr << "Caught Exception reading XML: " << e << std::endl;
            QDP_abort(1);
        }
    }

    // Function call
    void
    Inlinetests::operator()(unsigned long update_no,
                            XMLWriter &xml_out)
    {
        // If xml file not empty, then use alternate
        if (params.xml_file !="")
        {
            std::string xml_file = makeXMLFileName(params.xml_file, update_no);

            push(xml_out, "tests");
            write(xml_out, "update_no", update_no);
            write(xml_out, "xml_file", xml_file);
            pop(xml_out);

            XMLFileWriter xml(xml_file);
            func(update_no, xml);
        }
        else
        {
            func(update_no, xml_out);
        }
    }

    // Real work done here
    void
    Inlinetests::func(unsigned long update_no,
                      XMLWriter &xml_out)
    {
        START_CODE();

        QDPIO::cout << InlineTestsEnv::name << ": tests measurements" << std::endl;

        QDP::StopWatch snoop;
        snoop.reset();
        snoop.start();

        write(xml_out, "update_no", update_no);

        // Grab the gauge field
        multi1d<LatticeColorMatrix> u =
            TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix>>(params.named_obj.gauge_id);
        // Calculate some gauge invariant observables
        MesPlq(xml_out, "Observables", u);

        SftMom phases(0, true, 3);
        int tlength = phases.numSubsets();
        int xlength = Layout::lattSize()[2];
        int maxlength;
        multi2d<DComplex> hsum(4, tlength);
        multi1d<DComplex> wloop1(tlength);
        LatticeComplex tmp;
        multi1d<LatticeColorMatrix> ulink(4);
        LatticeColorMatrix u_tmp;

        for (int i = 0; i < 4; i++)
        {
            ulink[i] = u[i];
            int istep = 0;
            maxlength=xlength;
            if (i == 3)maxlength= tlength;
            for (int it = 0; it < maxlength; it++)
            {
                tmp = trace(ulink[i]);
                hsum[i][it] = sum(tmp);
                u_tmp = shift(ulink[i], FORWARD, i);
                ulink[i] = u[i] * u_tmp;
            }
        }

        push(xml_out, "tests"); // tests tag
        for (int i = 0; i < 4; i++)
        {
            write(xml_out, "wlink_v", hsum[i]);
        }
        pop(xml_out); // pop out tests tag

        snoop.stop();
        QDPIO::cout << InlineTestsEnv::name << ": total time = "
                    << snoop.getTimeInSeconds()
                    << " secs" << std::endl;

        QDPIO::cout << InlineTestsEnv::name << ": ran successfully" << std::endl;

        END_CODE();
    }
} // namespace Chroma