# Pumping_Test_Interpretation_with_Python

![Preview](https://numericalenvironmental.files.wordpress.com/2018/04/pumping_test_output.png)

This is a Python (version 2 or 3) script to aid in the interpretation of aquifer pumping tests using the Theis (1935) and Hantush and Jacob (1955) solutions for confined and leaky confined aquifers, respectively, as well as a simple numerical method-of-lines solution for an unconfined aquifer with wellbore storage. It is intended as a simple means for quick evaluation purposes and is not designed to be a comprehensive alternative to commercial codes such as AQTESOLVE [http://www.aqtesolv.com/]. The current version of the script does not perform automatic curve fitting but instead relies on the user to try out different combinations are aquifer/aquitard parameters and note the changes in the resulting model fits to the data. Additional discussion is provided in my blog post [https://numericalenvironmental.wordpress.com/2018/04/07/a-simple-set-of-analytical-and-numerical-solutions-in-python-for-pumping-test-interpretation/].
The following Python libraries are required:
* matplotlib (to plot the data and modeled drawdown curves)
* PyQt v5 (for the user interface)
* SciPy (for integration tools, including an ODE solver, as well as for special functions)

In addition, the following files must exist in the same folder as the script:
* aquifer.txt, a text file containing default aquifer parameter values
* well.txt, a text file describing pumping well characteristics
* transducer.txt, a text file with the monitor times and drawdown values representing the actual pumping test data
* pumping_test_interface.ui, and XML file associated with the user interface

To run the script, the user simply presses the “evaluate” button to generate a time-series drawdown plot comparing the data with the selected model output. The “update” button is pressed whenever a parameter value is changed on the user interface to permit the pumping test models to operate on the new parameter set. The “save” button updates the aquifer and well parameter text files. The following pumping test models are available, with the user checking all boxes that may apply:
* Confined (Theis): the classic Theis (1935) solution for confined aquifers
* Confined (wellbore storage): a numerical confined flow solution, including wellbore storage
* Leaky (Hantush and Jacob): the Hantush and Jacob (1955) leaky aquifer model, neglecting aquitard storage
* Unconfined (Theis): application of the Theis model to an unconfined aquifer, with specific yield used in place of specific storage. The impact of diminished saturated aquifer thickness near the wellbore is not considered.
* Unconfined (Dupuit): a numerical unconfined aquifer model which includes both wellbore storage and the impact of reduced saturated thickness on (horizontal) flow toward the well. This solution is for relatively thin aquifers that are fully intercepted by the well screen interval and does not address delayed yield. Note that excessive pumping with this model will generate a “dry” solution that will cause the script to crash. 

I'd appreciate hearing back from you if you find the code useful. Questions or comments are welcome at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

