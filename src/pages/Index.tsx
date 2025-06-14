
import React, { useState } from 'react';
import { Loader2, Search, FileText, Target, AlertCircle } from 'lucide-react';
import { Button } from '@/components/ui/button';
import { Textarea } from '@/components/ui/textarea';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Alert, AlertDescription } from '@/components/ui/alert';
import { Badge } from '@/components/ui/badge';
import { Separator } from '@/components/ui/separator';

const Index = () => {
  const [query, setQuery] = useState('');
  const [loading, setLoading] = useState(false);
  const [results, setResults] = useState<{
    summary: string;
    gaps: string[];
    status: string;
    message: string;
  } | null>(null);
  const [error, setError] = useState<string | null>(null);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    if (!query.trim()) return;

    setLoading(true);
    setError(null);
    setResults(null);

    try {
      // Replace with your actual Flask API endpoint
      const response = await fetch('/api/research', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ query: query.trim() }),
      });

      const data = await response.json();

      if (data.status === 'success') {
        setResults(data);
      } else {
        setError(data.message || 'An error occurred during analysis');
      }
    } catch (err) {
      setError('Failed to connect to the research service. Please try again.');
      console.error('API Error:', err);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="min-h-screen bg-gradient-to-br from-slate-50 to-blue-50">
      {/* Header */}
      <header className="bg-white shadow-sm border-b">
        <div className="max-w-6xl mx-auto px-4 py-6">
          <div className="flex items-center space-x-3">
            <div className="p-2 bg-blue-100 rounded-lg">
              <FileText className="h-6 w-6 text-blue-600" />
            </div>
            <div>
              <h1 className="text-2xl font-bold text-gray-900">
                NCD Literature Synthesis & Gap Identifier
              </h1>
              <p className="text-sm text-gray-600 mt-1">
                Intelligent analysis of non-communicable disease research literature
              </p>
            </div>
          </div>
        </div>
      </header>

      <main className="max-w-6xl mx-auto px-4 py-8">
        {/* Input Section */}
        <Card className="mb-8">
          <CardHeader>
            <CardTitle className="flex items-center space-x-2">
              <Search className="h-5 w-5" />
              <span>Research Query</span>
            </CardTitle>
          </CardHeader>
          <CardContent>
            <form onSubmit={handleSubmit} className="space-y-4">
              <div>
                <Textarea
                  placeholder="Enter your research query here... 

Examples:
• Community-based interventions for hypertension in rural Kenya
• Prevalence of type 2 diabetes among youth in Nairobi
• Cost-effectiveness of NCD prevention programs in sub-Saharan Africa"
                  value={query}
                  onChange={(e) => setQuery(e.target.value)}
                  className="min-h-[120px] text-base"
                  disabled={loading}
                />
              </div>
              <Button 
                type="submit" 
                size="lg" 
                disabled={!query.trim() || loading}
                className="w-full sm:w-auto"
              >
                {loading ? (
                  <>
                    <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                    Analyzing Literature...
                  </>
                ) : (
                  <>
                    <Search className="mr-2 h-4 w-4" />
                    Start Research Analysis
                  </>
                )}
              </Button>
            </form>
          </CardContent>
        </Card>

        {/* Error Display */}
        {error && (
          <Alert variant="destructive" className="mb-8">
            <AlertCircle className="h-4 w-4" />
            <AlertDescription>{error}</AlertDescription>
          </Alert>
        )}

        {/* Results Section */}
        {results && (
          <div className="space-y-6">
            {/* Success Message */}
            <Alert className="border-green-200 bg-green-50">
              <FileText className="h-4 w-4 text-green-600" />
              <AlertDescription className="text-green-800">
                {results.message}
              </AlertDescription>
            </Alert>

            {/* Literature Summary */}
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center space-x-2">
                  <FileText className="h-5 w-5 text-blue-600" />
                  <span>Literature Synthesis</span>
                </CardTitle>
              </CardHeader>
              <CardContent>
                <div className="prose prose-gray max-w-none">
                  <div className="bg-gray-50 rounded-lg p-6 text-gray-800 leading-relaxed whitespace-pre-wrap">
                    {results.summary}
                  </div>
                </div>
              </CardContent>
            </Card>

            {/* Research Gaps */}
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center space-x-2">
                  <Target className="h-5 w-5 text-orange-600" />
                  <span>Identified Research Gaps</span>
                  <Badge variant="secondary" className="ml-2">
                    {results.gaps.length} gaps found
                  </Badge>
                </CardTitle>
              </CardHeader>
              <CardContent>
                {results.gaps.length > 0 ? (
                  <div className="space-y-3">
                    {results.gaps.map((gap, index) => (
                      <div key={index} className="flex items-start space-x-3">
                        <div className="flex-shrink-0 w-6 h-6 bg-orange-100 rounded-full flex items-center justify-center text-sm font-medium text-orange-600 mt-0.5">
                          {index + 1}
                        </div>
                        <div className="flex-1 bg-orange-50 rounded-lg p-4 text-gray-800">
                          {gap}
                        </div>
                      </div>
                    ))}
                  </div>
                ) : (
                  <div className="text-center py-8 text-gray-500">
                    <Target className="h-12 w-12 mx-auto mb-3 text-gray-300" />
                    <p>No specific research gaps identified in the current analysis.</p>
                  </div>
                )}
              </CardContent>
            </Card>

            <Separator className="my-8" />

            {/* Query Information */}
            <Card className="bg-slate-50">
              <CardHeader>
                <CardTitle className="text-sm font-medium text-gray-600">
                  Analysis Details
                </CardTitle>
              </CardHeader>
              <CardContent>
                <div className="space-y-2">
                  <div>
                    <span className="text-sm font-medium text-gray-500">Query:</span>
                    <p className="text-sm text-gray-700 mt-1 bg-white rounded p-2">
                      {query}
                    </p>
                  </div>
                  <div>
                    <span className="text-sm font-medium text-gray-500">Status:</span>
                    <Badge variant="outline" className="ml-2">
                      {results.status}
                    </Badge>
                  </div>
                </div>
              </CardContent>
            </Card>
          </div>
        )}

        {/* Loading State */}
        {loading && (
          <Card>
            <CardContent className="py-12">
              <div className="text-center">
                <Loader2 className="h-8 w-8 animate-spin mx-auto mb-4 text-blue-600" />
                <h3 className="text-lg font-medium text-gray-900 mb-2">
                  Processing Your Research Query
                </h3>
                <p className="text-gray-600 max-w-md mx-auto">
                  Our AI agent is searching and analyzing relevant NCD literature. 
                  This may take a few moments...
                </p>
              </div>
            </CardContent>
          </Card>
        )}
      </main>

      {/* Footer */}
      <footer className="mt-16 py-8 border-t bg-white">
        <div className="max-w-6xl mx-auto px-4 text-center text-sm text-gray-500">
          <p>NCD Literature Synthesis & Gap Identifier - Powered by AI Research Assistant</p>
        </div>
      </footer>
    </div>
  );
};

export default Index;
