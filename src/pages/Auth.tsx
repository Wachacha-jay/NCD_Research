
import React, { useState } from 'react';
import { AuthForm } from '@/components/auth/AuthForm';
import { FileText } from 'lucide-react';

const Auth = () => {
  const [mode, setMode] = useState<'signin' | 'signup'>('signin');

  const toggleMode = () => {
    setMode(mode === 'signin' ? 'signup' : 'signin');
  };

  return (
    <div className="min-h-screen bg-gradient-to-br from-slate-50 to-blue-50 flex items-center justify-center p-4">
      <div className="w-full max-w-md">
        {/* Header */}
        <div className="text-center mb-8">
          <div className="flex items-center justify-center space-x-3 mb-4">
            <div className="p-2 bg-blue-100 rounded-lg">
              <FileText className="h-8 w-8 text-blue-600" />
            </div>
            <h1 className="text-2xl font-bold text-gray-900">
              NCD Research
            </h1>
          </div>
          <p className="text-gray-600">
            Access the Literature Synthesis & Gap Identifier
          </p>
        </div>

        <AuthForm mode={mode} onToggleMode={toggleMode} />
      </div>
    </div>
  );
};

export default Auth;
